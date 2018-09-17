import logging
import os

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.db.models import Q
from django.http import JsonResponse
from django.shortcuts import render, redirect
from django.utils import timezone
from django.utils.safestring import mark_safe
from sendfile import sendfile

from biostar.accounts.models import Profile, User
from biostar.forum import views as forum_views
from biostar.forum.models import Post
from biostar.utils.shortcuts import reverse
from biostar.utils.decorators import ajax_error_wrapper

from . import tasks, auth, forms, util, const
from .decorators import object_access
from .diffs import color_diffs
from .models import (Project, Data, Analysis, Job, Access)

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


def index(request):
    context = dict()
    return render(request, 'index.html', context)


@user_passes_test(lambda u: u.is_superuser)
def site_admin(request):
    '''
    Administrative view. Lists the admin project and job.
    '''
    projects = Project.objects.all()
    context = dict(projects=projects)
    return render(request, 'admin_index.html', context=context)


@login_required
def recycle_bin(request):
    "Recycle bin view for a user"

    # Only searches projects you have access.
    all_projects = auth.get_project_list(user=request.user)

    del_data = Data.objects.get_deleted(project__in=all_projects,
                                        owner=request.user).order_by("date")

    del_recipes = Analysis.objects.get_deleted(project__in=all_projects,
                                               owner=request.user).order_by("date")

    del_jobs = Job.objects.get_deleted(project__in=all_projects,
                                       owner=request.user).order_by("date")

    context = dict(jobs=del_jobs, data=del_data, recipes=del_recipes)

    return render(request, 'recycle_bin.html', context=context)


@login_required
def recipe_mod(request):
    "Shows all recipes that are under review."

    user = request.user

    if not user.profile.is_manager:
        msg = mark_safe("You have to be a <b>manager</b> to view this page.")
        messages.error(request, msg)
        return redirect("/")

    context = dict()
    return render(request, 'recipe_mod.html', context)


def clear_clipboard(request, uid):
    "Clear copied objects held in clipboard."

    next_url = request.GET.get("next", reverse("project_view", kwargs=dict(uid=uid)))
    board = request.GET.get("board")

    if board:
        request.session[board] = []

    return redirect(next_url)


def get_access(request, project):
    user = request.user if request.user.is_authenticated else None
    user_access = Access.objects.filter(project=project, user=user).first()
    # Current users access
    user_access = user_access or Access(access=Access.NO_ACCESS)

    # Users already with access to current project
    user_list = [a.user for a in project.access_set.all() if a.access > Access.NO_ACCESS]

    return user_access, user_list


@object_access(type=Project, access=Access.OWNER_ACCESS, url='data_list')
def project_users(request, uid):
    """
    Manage project users
    """
    project = Project.objects.filter(uid=uid).first()
    user_access, user_list = get_access(request, project)
    label = lambda x: f"<span class='ui green tiny label'>{x}</span>"

    # Search query separate for users.
    q = request.GET.get("q", "")
    form = forms.ChangeUserAccess()

    if request.method == "POST":
        form = forms.ChangeUserAccess(data=request.POST)

        # User needs to be authenticated and have admin access to make any changes.
        if form.is_valid() and request.user.is_authenticated:
            user, access = form.save()
            msg = f"Changed <b>{user.first_name}</b>'s access to {label(access.get_access_display())}"
            messages.success(request, mark_safe(msg))
            return redirect(reverse("project_users", request=request, kwargs=dict(uid=project.uid)))

    # Users that have been searched for.
    targets = User.objects.filter(Q(email__contains=q) | Q(first_name__contains=q)) if q else []
    current = forms.access_forms(users=user_list, project=project, exclude=[request.user])
    results = forms.access_forms(users=targets, project=project, exclude=[request.user])
    context = dict(current=current, project=project, results=results, form=form, activate='User Management',
                   q=q, user_access=user_access)
    counts = get_counts(project)
    context.update(counts)
    return render(request, "project_users.html", context=context)


def project_list(request):
    projects = auth.get_project_list(user=request.user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")

    context = dict(projects=projects)
    return render(request, "project_list.html", context)


@object_access(type=Project, access=Access.READ_ACCESS)
def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """

    return project_view(request=request, uid=uid, template_name="data_list.html",
                        active='data', show_summary=True)


@object_access(type=Project, access=Access.READ_ACCESS)
def discussion_list(request, uid):
    posts = Post.objects.get_discussions(project__uid=uid, type__in=Post.TOP_LEVEL).order_by("-pk").all()

    context = dict(posts=posts)
    return project_view(request=request, uid=uid, template_name="discussion_list.html",
                        active='discussion', extra_context=context)


@object_access(type=Post, access=Access.READ_ACCESS)
def discussion_subs(request, uid):
    next_url = reverse("discussion_view", request=request, kwargs=dict(uid=uid))
    return forum_views.subs_action(request=request, uid=uid, next=next_url)


@object_access(type=Project, access=Access.WRITE_ACCESS, login_required=True, url="discussion_list")
def discussion_create(request, uid):
    project = Project.objects.filter(uid=uid).first()

    template = "discussion_create.html"

    context = dict(project=project, activate='Start a Discussion')
    counts = get_counts(project)
    context.update(counts)

    allowed_posts = [
        Post.QUESTION, Post.COMMENT, Post.ANSWER, Post.TOOL,
        Post.TUTORIAL, Post.DATA, Post.NEWS
    ]
    filter_function = lambda x: x[0] in allowed_posts

    return forum_views.post_create(request=request, template=template, extra_context=context,
                                   url="discussion_view", project=project, filter_func=filter_function)


@object_access(type=Post, access=Access.READ_ACCESS)
def discussion_view(request, uid):
    template = "discussion_view.html"
    # Get the parents info
    obj = Post.objects.get_discussions(uid=uid).first()

    project = obj.root.project
    comment_url = reverse("discussion_comment")

    context = dict(project=project, activate="Discussion", comment_url=comment_url)
    counts = get_counts(project)
    context.update(counts)

    return forum_views.post_view(request=request, template=template, extra_context=context,
                                 url="discussion_view", uid=uid)


@object_access(type=Project, access=Access.READ_ACCESS)
def recipe_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """

    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes')


def job_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


def get_counts(project):
    data_count = project.data_set.count()
    recipe_count = project.analysis_set.count()
    result_count = project.job_set.count()
    discussion_count = Post.objects.get_discussions(project=project,
                                                    type__in=Post.TOP_LEVEL).count()

    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count,
        discussion_count=discussion_count
    )


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, uid, template_name="recipe_list.html", active='recipes', show_summary=None,
                 extra_context={}):
    project = Project.objects.filter(uid=uid).first()
    # Show counts for the project.
    counts = get_counts(project)

    # Select all the data in the project.
    data_list = project.data_set.order_by("-sticky", "-date").all()
    recipe_list = project.analysis_set.order_by("-sticky", "-date").all()
    job_list = project.job_set.order_by("-sticky", "-date").all()

    # Filter job results by analysis
    recipe_filter = request.GET.get('filter', '')
    if recipe_filter:
        job_list = job_list.filter(analysis__uid=recipe_filter)

    context = dict(project=project, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, recipe_filter=recipe_filter, show_summary=show_summary)
    context.update(counts)
    context.update(extra_context)

    return render(request, template_name, context)


@object_access(type=Project, access=Access.OWNER_ACCESS, url='data_list')
def project_edit(request, uid):
    "Edit meta-data associated with a project."

    project = Project.objects.filter(uid=uid).first()
    form = forms.ProjectForm(instance=project)
    if request.method == "POST":
        form = forms.ProjectForm(request.POST, request.FILES, instance=project)
        if form.is_valid():
            form.save()
            return redirect(reverse("project_view", request=request, kwargs=dict(uid=project.uid)))

    context = dict(project=project, form=form)
    return render(request, "project_edit.html", context=context)


@login_required
def project_create(request):
    """
    View used create an empty project belonging to request.user.
    Input is validated with a form and actual creation is routed through auth.create_project.
    """
    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = forms.ProjectForm(initial=initial)

    if request.method == "POST":
        # create new projects here ( just populates metadata ).
        form = forms.ProjectForm(request.POST, request.FILES)
        if form.is_valid():
            project = form.custom_save(owner=request.user)
            return redirect(reverse("project_view", request=request, kwargs=dict(uid=project.uid)))

    context = dict(form=form)
    return render(request, "project_create.html", context=context)


def ajax_copy(request, modeltype, msg="Copied!", board=None):

    user = request.user
    data_uid = request.GET.get("data_uid")
    instance = modeltype.objects.get_all(uid=data_uid).first()

    if instance is None:
        entry = Access(access=Access.NO_ACCESS)
    else:
        entry = Access.objects.filter(user=user, project=instance.project).first()

    if entry.access >= Access.READ_ACCESS or instance.project.is_public:
        current = request.session.get(board, [])
        current.append(instance.uid)
        # No duplicates in clipboard
        current = list(set(current))
        request.session[board] = current
        status = "success"
        msg += f" Clipboard has {len(current)} object(s)."
    else:
        msg = status = "error"
        request.session[board] = []

    response = JsonResponse({"msg": msg, "status": status})
    return response


@ajax_error_wrapper(method="GET")
def ajax_job_copy(request):

    msg = "Copied to results clipboard!"
    return ajax_copy(modeltype=Job, request=request, msg=msg, board=const.RESULTS_CLIPBOARD)


@ajax_error_wrapper(method="GET")
def ajax_data_copy(request):

    msg = "Copied to data clipboard!"
    return ajax_copy(modeltype=Data, request=request, msg=msg, board=const.DATA_CLIPBOARD)


@ajax_error_wrapper(method="GET")
def ajax_recipe_copy(request):

    msg = "Copied to recipe clipboard!"
    return ajax_copy(modeltype=Analysis, request=request, msg=msg, board=const.RECIPE_CLIPBOARD)


@object_access(type=Project, access=Access.WRITE_ACCESS, url="recipe_list")
def recipe_paste(request, uid):
    """Used to paste recipes in a clipboard as a new recipes."""

    project = Project.objects.filter(uid=uid).first()

    clipboard = request.session.get(const.RECIPE_CLIPBOARD, [])

    for uid in clipboard:
        instance = Analysis.objects.get_all(uid=uid).first()

        if instance:
            name, summary, text = instance.name, instance.summary, instance.text
            new_recipe = auth.create_analysis(project=project, json_text=instance.json_text,
                                              template=instance.template, summary=summary, user=request.user, name=name,
                                              text=text, stream=instance.image, security=instance.security)
            # Ensure the diff gets inherited.
            new_recipe.last_valid = instance.last_valid
            new_recipe.save()

    request.session[const.RECIPE_CLIPBOARD] = []
    messages.success(request, mark_safe(f"Pasted <b>{len(clipboard)} recipes</b>  in clipboard"))
    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.WRITE_ACCESS, url="data_list")
def data_paste(request, uid):
    """Used to paste objects in results and data clipboards as a Data object."""

    project = Project.objects.filter(uid=uid).first()
    owner = request.user
    board = request.GET.get("board")
    clipboard = request.session.get(board, [])

    for datauid in clipboard:

        if board == const.DATA_CLIPBOARD:
            obj = Data.objects.get_all(uid=datauid).first()
            dtype = obj.type

        else:
            obj = Job.objects.get_all(uid=datauid).first()
            dtype = "DATA"

        if obj:
            paths = [n.path for n in os.scandir(obj.get_data_dir())]
            auth.create_data(project=project, paths=paths, user=owner,
                             name=obj.name, type=dtype, summary=obj.summary)

    request.session[board] = []
    messages.success(request, "Pasted data in clipboard")
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Data, access=Access.READ_ACCESS)
def data_view(request, uid):
    "Show information specific to each data."

    data = Data.objects.get_all(uid=uid).first()
    project = data.project

    context = dict(data=data, project=project, activate='Selected Data')
    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.OWNER_ACCESS, url='data_view')
def data_edit(request, uid):
    """
    Edit meta-data associated with Data.
    """

    data = Data.objects.get_all(uid=uid).first()
    form = forms.DataEditForm(instance=data, initial=dict(type=data.type), user=request.user)

    if request.method == "POST":
        form = forms.DataEditForm(data=request.POST, instance=data, user=request.user, files=request.FILES)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", request=request, kwargs=dict(uid=data.uid)))
        print(form.errors)
    context = dict(data=data, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.WRITE_ACCESS, url='data_list')
def data_upload(request, uid):
    "Data upload view routed through auth.create_data."

    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    form = forms.DataUploadForm(user=owner, project=project)

    if request.method == "POST":
        form = forms.DataUploadForm(data=request.POST, files=request.FILES, user=owner, project=project)

        if form.is_valid():
            data = form.save()
            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", request=request, kwargs={'uid': project.uid}))

    context = dict(project=project, form=form, activate="Add Data")

    counts = get_counts(project)

    context.update(counts)

    return render(request, 'data_upload.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER)
def recipe_view(request, uid):
    """
    Returns an analysis view based on its id.
    """
    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project
    context = dict(recipe=recipe, project=project, activate='View Recipe')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "recipe_view.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view')
def recipe_code_view(request, uid):
    """
    Returns an analysis code view based on its id.
    """
    user = request.user
    recipe = Analysis.objects.get_all(uid=uid).first()
    form = forms.RecipeCodeEdit(user=user, recipe=recipe)

    if request.method == "POST":
        form = forms.RecipeCodeEdit(data=request.POST, user=user, recipe=recipe)
        if form.is_valid():

            template = form.cleaned_data.get('template').strip()

            if template != recipe.template:
                recipe.template = template
                recipe.security = auth.authorize_analysis(user=user, recipe=recipe)
                recipe.diff_author = user
                recipe.diff_date = timezone.now()
                recipe.save()
                messages.success(request, f"The recipe has been updated.")
            else:
                messages.info(request, f"The recipe has not been modified.")

            return redirect(reverse("recipe_code_view", request=request, kwargs={'uid': recipe.uid}))

    project = recipe.project
    context = dict(recipe=recipe, project=project, activate='Recipe Code', form=form)

    counts = get_counts(project)
    context.update(counts)

    return render(request, "recipe_code_view.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view', show_deleted=False)
def recipe_run(request, uid):
    """
    View used to execute recipes and start a 'Queued' job.
    """

    analysis = Analysis.objects.filter(uid=uid).first()
    project = analysis.project

    if request.method == "POST":
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data,
                                     data=request.POST)

        if form.is_valid():

            # The desired name of for the results.
            name = form.cleaned_data.get("name")

            # Generates the JSON data from the bound form field.
            json_data = form.fill_json_data()

            # Create the job from the json.
            state = Job.SPOOLED if tasks.HAS_UWSGI else Job.QUEUED
            job = auth.create_job(analysis=analysis, user=request.user, json_data=json_data, name=name, state=state)

            # Spool the job right away if UWSGI exists.
            if tasks.HAS_UWSGI:
                tasks.execute_job.spool(job_id=job.id)

            return redirect(reverse("job_list", request=request, kwargs=dict(uid=project.uid)))
    else:
        initial = dict(name=f"Results for: {analysis.name}")
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, form=form, activate='Run Recipe')
    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER, url='recipe_view')
def recipe_code_edit(request, uid):
    """
    Displays and allows edit on a recipe code.

    Since we allow a preview for un-authenticated users thus the view
    is more complicated than a typical DJANGO form handler.
    """
    user = request.user

    # There has to be a recipe to work with.
    analysis = Analysis.objects.get_all(uid=uid).first()
    project = analysis.project
    name = analysis.name

    if request.method == "POST":
        form = forms.EditCode(user=user, project=project, data=request.POST)

        if form.is_valid():

            # Templates.
            template = form.cleaned_data['template']

            # Preview action will let the form cascade through.
            save = form.cleaned_data['action'] == 'SAVE'

            analysis.json_text = form.cleaned_data['json']

            # Changes to template will require a review ( only when saving ).
            if auth.template_changed(analysis=analysis, template=template) and save:
                analysis.security = Analysis.UNDER_REVIEW
                analysis.diff_author = user
                analysis.diff_date = timezone.now()

            # Staff members will automatically get authorized.
            if user.is_staff:
                analysis.security = Analysis.AUTHORIZED

            # Set the new template.
            analysis.template = template

            # Only the SAVE action commits the changes on the analysis.
            if save:
                analysis.save()
                messages.info(request, "The recipe has been updated.")
                return redirect(reverse("recipe_view", request=request, kwargs=dict(uid=analysis.uid)))
    else:
        # This gets triggered on a GET request.
        initial = dict(template=analysis.template, json=analysis.json_text)
        form = forms.EditCode(user=user, project=project, initial=initial)

    # Bind the JSON to the form.
    recipe = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data,
                                   initial=dict(name=name))

    # This generates a "fake" unsaved job.
    # Needs to fill in a few runtime only settings.
    mock_json = analysis.json_data
    for key, value in mock_json.items():
        value['toc'] = f"{key}-filelist.txt"
    job = auth.create_job(analysis=analysis, json_data=mock_json, save=False)

    # Create the script for the "fake" job.
    data, script = auth.generate_script(job)

    # Populate the context.
    context = dict(project=project, analysis=analysis, form=form, script=script, recipe=recipe)
    return render(request, 'recipe_edit_code.html', context)


@object_access(type=Project, access=Access.WRITE_ACCESS, url='recipe_list')
def recipe_create(request, uid):
    """
    Create recipe with empty template and json spec
    """

    project = Project.objects.filter(uid=uid).first()
    form = forms.RecipeForm(initial=dict(name="New Recipe"))

    if request.method == "POST":
        form = forms.RecipeForm(data=request.POST, files=request.FILES)

        if form.is_valid():
            # Recipe is authorized since the template is empty at this point.
            security = Analysis.AUTHORIZED
            name = form.cleaned_data["name"]
            text = form.cleaned_data["text"]
            summary = form.cleaned_data["summary"]
            stream = form.cleaned_data["image"]
            sticky = form.cleaned_data["sticky"]
            uid = form.cleaned_data["uid"]

            recipe = auth.create_analysis(project=project, json_text="{}", template="",
                                          user=request.user, summary=summary, name=name, text=text,
                                          security=security, stream=stream, sticky=sticky, uid=uid)
            recipe.save()
            messages.success(request, "Recipe created")

            return redirect(reverse('recipe_list', request=request, kwargs=dict(uid=project.uid)))
    # The url to submit to.
    action_url = reverse('recipe_create', request=request, kwargs=dict(uid=project.uid))
    context = dict(project=project, form=form, action_url=action_url, name="New Recipe")

    return render(request, 'recipe_edit.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER, url='recipe_view')
def recipe_diff(request, uid):
    """
    View used to show diff in template and authorize it.
    """
    recipe = Analysis.objects.get_all(uid=uid).first()
    differ = auth.template_changed(template=recipe.last_valid, analysis=recipe)
    differ = color_diffs(differ)

    form = forms.RecipeDiff(recipe=recipe, request=request, user=request.user)

    if request.method == "POST":
        form = forms.RecipeDiff(recipe=recipe, user=request.user, data=request.POST,
                                request=request)
        if form.is_valid():
            form.save()
            return redirect(reverse('recipe_view', request=request, kwargs=dict(uid=recipe.uid)))

    context = dict(activate="Recent Template Change", project=recipe.project, recipe=recipe,
                   diff=mark_safe(''.join(differ)), form=form)

    counts = get_counts(recipe.project)
    context.update(counts)

    return render(request, "recipe_diff.html", context=context)


@object_access(type=Analysis, access=Access.OWNER_ACCESS, url='recipe_view')
def recipe_edit(request, uid):
    "Edit meta-data associated with a recipe."

    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project
    action_url = reverse('recipe_edit', request=request, kwargs=dict(uid=recipe.uid))
    form = forms.RecipeForm(instance=recipe)

    if request.method == "POST":
        form = forms.RecipeForm(data=request.POST, files=request.FILES, instance=recipe)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", request=request, kwargs=dict(uid=recipe.uid)))

    context = dict(analysis=recipe, project=project, form=form, action_url=action_url,
                   name=recipe.name)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Job, access=Access.OWNER_ACCESS, url="job_view")
def job_edit(request, uid):
    "Edit meta-data associated with a job."

    job = Job.objects.get_all(uid=uid).first()
    project = job.project
    form = forms.JobEditForm(instance=job)

    if request.method == "POST":
        form = forms.JobEditForm(data=request.POST, files=request.FILES, instance=job)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", request=request, kwargs=dict(uid=job.uid)))

    context = dict(job=job, project=project, form=form)
    return render(request, 'job_edit.html', context)


@object_access(type=Analysis, access=Access.OWNER_ACCESS, url="recipe_view")
def recipe_delete(request, uid):

    recipe = Analysis.objects.get_all(uid=uid).first()

    auth.delete_object(obj=recipe, request=request)

    return redirect(reverse("recipe_list", kwargs=dict(uid=recipe.project.uid)))


@object_access(type=Job, access=Access.OWNER_ACCESS, url="job_view")
def job_delete(request, uid):

    job = Job.objects.get_all(uid=uid).first()

    running_job = job.state == Job.RUNNING and not job.deleted

    if running_job:
        messages.error(request, "Can not delete a running job. Wait until it finishes.")
        return redirect(job.url())

    auth.delete_object(obj=job, request=request)
    return redirect(reverse("job_list", kwargs=dict(uid=job.project.uid)))


@object_access(type=Data, access=Access.OWNER_ACCESS, url="data_view")
def data_delete(request, uid):

    data = Data.objects.get_all(uid=uid).first()

    auth.delete_object(obj=data, request=request)

    return redirect(reverse("data_list", kwargs=dict(uid=data.project.uid)))


@object_access(type=Job, access=Access.READ_ACCESS)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.get_all(uid=uid).first()
    project = job.project

    # The path is a GET parameter
    path = request.GET.get('path', "")

    # The job rooth directory
    root = job.path

    # Get the target directory.
    abspath = join(job.path, path)

    # Generate the files
    try:
        files = util.scan_files(abspath=abspath, relpath=path, root=root)
    except Exception as exc:
        messages.error(request, f"{exc}")
        files = []

    context = dict(job=job, project=project, activate='View Result', files=files, path=path)

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def file_serve(request, path, obj):
    """
    Authenticates access through decorator before serving file.
    """

    # Get the object that corresponds to the entry.
    root = obj.get_data_dir()

    # This will turn into an absolute path.
    file_path = join(root, path)

    # Ensure only files in the object root can be accessed.
    if not file_path.startswith(root) or not os.path.isfile(file_path):
        msg = "Invalid path." if (not file_path.startswith(root)) else f"File not found: {path}"
        messages.error(request, msg)
        return redirect(obj.url())

    # The response will be the file content.
    mimetype = auth.guess_mimetype(fname=path)

    # Get the filesize in Mb
    size = os.path.getsize(file_path) / 1024 / 1024

    # This behavior can be further customized in front end webserver.
    if size < 20:
        # Return small files in the browser if possible.
        data = sendfile(request, file_path, mimetype=mimetype)
    else:
        # Trigger a file download for bigger files.
        fname = os.path.basename(file_path)
        data = sendfile(request, file_path, attachment=True, attachment_filename=fname, mimetype=mimetype)

    return data


@object_access(type=Data, access=Access.READ_ACCESS, url='data_navigate')
def data_serve(request, uid, path):
    """
    Serves files from a data directory.
    """
    obj = Data.objects.get_all(uid=uid).first()
    return file_serve(request=request, path=path, obj=obj)


def job_serve(request, uid, path):
    """
    Serves files from a job directory.
    """
    obj = Job.objects.get_all(uid=uid).first()

    if obj:
        return file_serve(request=request, path=path, obj=obj)
    else:
        messages.error(request, "Object does not exist")
        return redirect("/")
