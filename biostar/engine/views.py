import logging
import os

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.db.models import Sum
from django.shortcuts import render, redirect
from django.template import Template, Context
from django.db.models import Count
from django.utils import timezone
from django.utils.safestring import mark_safe
from ratelimit.decorators import ratelimit
from sendfile import sendfile
from django.http import HttpResponse
from django.conf import settings
from django.db.models import Q, Count

from biostar.accounts.models import User
from biostar.forum import views as forum_views
from biostar.forum.models import Post
from biostar.utils.shortcuts import reverse
from . import tasks, auth, forms, const, util, search
from .decorators import read_access, write_access
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
    user = request.user

    if user.is_superuser:
        # Super users get access to all deleted objects.
        projects = Project.objects.get_all()
        query_dict = dict(project__in=projects)
    else:
        # Only searches projects user have access.
        projects = auth.get_project_list(user=user, include_deleted=True)
        query_dict = dict(project__in=projects, owner=user)

    projects = projects.filter(deleted=True).order_by("date")
    projects = annotate_projects(projects)
    data = Data.objects.get_deleted(**query_dict).order_by("date")
    recipes = Analysis.objects.get_deleted(**query_dict).order_by("date")
    jobs = Job.objects.get_deleted(**query_dict).order_by("date")

    context = dict(jobs=jobs, data=data, recipes=recipes, projects=projects)

    return render(request, 'recycle_bin.html', context=context)


@write_access(type=Project, fallback_view="project_view")
def project_delete(request, uid):

    project = Project.objects.get_all(uid=uid).first()
    project.deleted = not project.deleted
    project.save()

    msg = f"Project:{project.name} successfully "
    msg += "deleted!" if project.deleted else "restored!"

    messages.success(request, msg)

    return redirect(reverse("project_list_private"))


def clear_clipboard(request, uid):
    "Clear copied objects held in clipboard."

    next_url = request.GET.get("next", reverse("project_view", kwargs=dict(uid=uid)))
    board = request.GET.get("board")
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})

    if clipboard.get(board):
        clipboard[board] = []
        request.session.update({settings.CLIPBOARD_NAME: clipboard})

    return redirect(next_url)


def search_bar(request):

    results = search.search(request=request)

    # Indicate to users that minimum character needs to be met.
    query_lenth = len(request.GET.get("q", "").strip())
    min_length = query_lenth > settings.SEARCH_CHAR_MIN

    # Indicate to users that there are no results for search.
    current_results = len([inner for outer in results.values() for inner in outer])
    no_results = min_length and current_results == 0

    context = dict(results=results, query=request.GET.get("q", "").strip(),
                   min_length=min_length, no_results=no_results)

    return render(request, "search.html", context)


def get_access(request, project):
    user = request.user if request.user.is_authenticated else None
    user_access = Access.objects.filter(project=project, user=user).first()
    # Current users access
    user_access = user_access or Access(access=Access.NO_ACCESS)

    # Users already with access to current project
    user_list = [a.user for a in project.access_set.all() if a.access > Access.NO_ACCESS]

    return user_access, user_list


@write_access(type=Project, fallback_view="data_list")
def project_users(request, uid):
    """
    Manage project users
    """
    project = Project.objects.get_all(uid=uid).first()
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


@read_access(type=Project)
def project_info(request, uid):

    user = request.user

    project = Project.objects.get_all(uid=uid).first()

    # Show counts for the project.
    counts = get_counts(project)

    # Who has write access
    write_access = auth.has_write_access(user=user, project=project)

    context = dict(project=project, active="info", write_access=write_access)
    context.update(counts)

    return render(request, "project_info.html", context)


def annotate_projects(projects):
    projects = projects.annotate(data_count=Count('data', distinct=True, filter=Q(deleted=False)),
                                 job_count=Count('job', distinct=True, filter=Q(deleted=False)),
                                 recipe_count=Count('analysis', distinct=True, filter=Q(deleted=False)))
    return projects


def project_list_private(request):
    """Only list private projects belonging to a user."""

    projects = auth.get_project_list(user=request.user, include_public=False)

    empty_msg = "No projects found."
    if request.user.is_anonymous:
        projects = []
        empty_msg = mark_safe(f"You need to <a href={reverse('login')}> log in</a> to view your projects.")
    else:
        projects = projects.order_by("rank", "-date", "-lastedit_date", "-id")
        projects = annotate_projects(projects)

    context = dict(projects=projects, msg=empty_msg)

    return render(request, "project_list_private.html", context)


def project_list_public(request):
    """Only list public projects."""

    projects = auth.get_project_list(user=request.user)
    # Exclude private projects
    projects = projects.exclude(privacy=Project.PRIVATE)
    projects = projects.order_by("rank", "-date", "-lastedit_date", "-id")
    projects = annotate_projects(projects)

    context = dict(projects=projects)

    return render(request, "project_list_public.html", context)


def project_list(request):

    if request.user.is_authenticated:
        # Return private projects when user is logged in.
        return project_list_private(request)
    else:
        return project_list_public(request)


@read_access(type=Project)
def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """

    return project_view(request=request, uid=uid, template_name="data_list.html",
                        active='data', show_summary=True)


@read_access(type=Project)
def discussion_list(request, uid):
    posts = Post.objects.get_discussions(project__uid=uid, type__in=Post.TOP_LEVEL).order_by("-pk").all()

    context = dict(posts=posts)
    return project_view(request=request, uid=uid, template_name="discussion_list.html",
                        active='discussion', extra_context=context)


@read_access(type=Post)
def discussion_subs(request, uid):
    next_url = reverse("discussion_view", request=request, kwargs=dict(uid=uid))

    return forum_views.subs_action(request=request, uid=uid, next=next_url)


@write_access(type=Project, fallback_view="discussion_list")
def discussion_create(request, uid):
    project = Project.objects.get_all(uid=uid).first()

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


@read_access(type=Post)
def discussion_view(request, uid):
    template = "discussion_view.html"
    # Get the parents info
    obj = Post.objects.get_discussions(uid=uid).first()

    project = obj.root.project
    sub_url = reverse("discussion_subs", kwargs=dict(uid=obj.uid))
    next_url = reverse("discussion_view", kwargs=dict(uid=obj.uid))

    context = dict(project=project, activate="Discussion", sub_url=sub_url, next_url=next_url)
    counts = get_counts(project)
    context.update(counts)

    return forum_views.post_view(request=request, template=template, extra_context=context,
                                 url="discussion_view", uid=uid)


@read_access(type=Project)
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
    discussion_count = 0

    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count,
        discussion_count=discussion_count
    )


@read_access(type=Project)
def project_view(request, uid, template_name="project_info.html", active='info', show_summary=None,
                 extra_context={}):
    """
    This view handles the project info, data list, recipe list, result list views.
    """

    # The user making the request
    user = request.user

    # The project that is viewed.
    project = Project.objects.get_all(uid=uid).first()

    # Select all the data in the project.
    data_list = project.data_set.order_by("rank", "-date").all()
    recipe_list = project.analysis_set.order_by("rank", "-date").all()

    # Annotate each recipe with the number of jobs it has.
    recipe_list = recipe_list.annotate(job_count=Count("job", filter=Q(job__deleted=False)))

    job_list = project.job_set.order_by("-date").all()

    # Filter job results by analysis
    filter_uid = request.GET.get('filter', '')
    recipe_filter = Analysis.objects.filter(uid=filter_uid).first()

    # The recipe filter exists
    if recipe_filter:
        job_list = job_list.filter(analysis=recipe_filter)

    # Add related content.
    job_list = job_list.select_related("analysis")

    # Who has write access
    write_access = auth.has_write_access(user=user, project=project)

    # Build the context for the project.
    context = dict(project=project, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, recipe_filter=recipe_filter, write_access=write_access)

    # Compute counts for the project.
    counts = get_counts(project)

    # Update conext with the counts.
    context.update(counts)

    # Add any extra context that may come from parameters.
    context.update(extra_context)

    return render(request, template_name, context)


@write_access(type=Project, fallback_view="data_list")
def project_edit(request, uid):
    "Edit meta-data associated with a project."

    project = Project.objects.get_all(uid=uid).first()
    form = forms.ProjectForm(instance=project, request=request)
    if request.method == "POST":
        form = forms.ProjectForm(data=request.POST, files=request.FILES, instance=project, request=request)
        if form.is_valid():
            project = form.save()
            Project.objects.get_all(uid=uid).update(lastedit_user=request.user)
            return redirect(reverse("project_view", request=request, kwargs=dict(uid=project.uid)))

    context = dict(project=project, form=form)
    return render(request, "project_edit.html", context=context)


@login_required
@ratelimit(key='ip', rate='5/h', block=True, method=ratelimit.UNSAFE)
def project_create(request):
    """
    View used create an empty project belonging to request.user.
    Input is validated with a form and actual creation is routed through auth.create_project.
    """
    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = forms.ProjectForm(initial=initial, request=request, create=True)

    if request.method == "POST":
        # create new projects here ( just populates metadata ).
        form = forms.ProjectForm(request=request, data=request.POST, create=True, files=request.FILES)
        if form.is_valid():
            project = form.custom_save(owner=request.user)
            return redirect(reverse("project_view", request=request, kwargs=dict(uid=project.uid)))

    context = dict(form=form)
    return render(request, "project_create.html", context=context)


@read_access(type=Data)
def data_copy(request, uid):
    data = Data.objects.get_all(uid=uid).first()
    next_url = request.GET.get("next", reverse("data_list", kwargs=dict(uid=data.project.uid)))

    auth.copy_uid(request=request, instance=data, board=const.DATA_CLIPBOARD)

    return redirect(next_url)


@read_access(type=Analysis)
def recipe_copy(request, uid):
    recipe = Analysis.objects.get_all(uid=uid).first()
    next_url = request.GET.get("next", reverse("recipe_list", kwargs=dict(uid=recipe.project.uid)))

    auth.copy_uid(request=request, instance=recipe, board=const.RECIPE_CLIPBOARD)

    return redirect(next_url)


@read_access(type=Job)
def job_copy(request, uid):
    job = Job.objects.get_all(uid=uid).first()
    next_url = request.GET.get("next", reverse("job_list", kwargs=dict(uid=job.project.uid)))

    auth.copy_uid(request=request, instance=job, board=const.RESULTS_CLIPBOARD)

    return redirect(next_url)


@read_access(type=Data)
def data_file_copy(request, uid, path):

    # Get the root data where the file exists
    data = Data.objects.get_all(uid=uid).first()
    fullpath = os.path.join(data.get_data_dir(), path)
    auth.copy_file(request=request, fullpath=fullpath)

    return redirect(reverse("data_view", kwargs=dict(uid=uid)))


@read_access(type=Job)
def job_file_copy(request, uid, path):

    # Get the root data where the file exists
    job = Job.objects.get_all(uid=uid).first()
    fullpath = os.path.join(job.get_data_dir(), path)

    auth.copy_file(request=request, fullpath=fullpath)

    return redirect(reverse("job_view", kwargs=dict(uid=uid)))


@write_access(type=Project, fallback_view="recipe_list")
def recipe_paste(request, uid):
    """
    Pastes recipes from clipboard as a new recipes.
    """

    # The user performing the action.
    user = request.user

    # The project the paste will use.
    project = Project.objects.get_all(uid=uid).first()

    # Contains the uids for the recipes that are to be copied.
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    recipe_uids = clipboard.get(const.RECIPE_CLIPBOARD, [])

    # Select valid recipe uids.
    recipes = [Analysis.objects.get_all(uid=uid).first() for uid in recipe_uids]

    # Keep existing recipes.
    recipes = filter(None, recipes)

    # The copy function for each recipe.
    def copy(instance):
        recipe = auth.create_analysis(project=project, user=user,
                                      json_text=instance.json_text,
                                      template=instance.template,
                                      name=instance.name, text=instance.text, stream=instance.image)
        return recipe

    # The list of new object created by the copy.
    new_recipes = list(map(copy, recipes))

    # Reset the session.
    clipboard[const.RECIPE_CLIPBOARD] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})

    # Notification after paste.
    messages.success(request, mark_safe(f"Pasted <b>{len(new_recipes)} recipes</b>  in clipboard"))

    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@write_access(type=Project, fallback_view="data_list")
def data_paste(request, uid):
    """Used to paste objects in results and data clipboards as a Data object."""
    project = Project.objects.get_all(uid=uid).first()
    owner = request.user
    board = request.GET.get("board")
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    data_clipboard = clipboard.get(board, [])

    for datauid in data_clipboard:

        if board == const.DATA_CLIPBOARD:
            obj = Data.objects.get_all(uid=datauid).first()
            dtype = obj.type
        else:
            obj = Job.objects.get_all(uid=datauid).first()
            dtype = "DATA"

        if obj:
            auth.create_data(project=project, path=obj.get_data_dir(), user=owner, name=obj.name,
                             type=dtype, text=obj.text)

    clipboard[board] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})
    messages.success(request, "Pasted data in clipboard")
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@write_access(type=Project, fallback_view="data_list")
def file_paste(request, uid):

    project = Project.objects.get_all(uid=uid).first()
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    file_clipboard = clipboard.get(const.FILES_CLIPBOARD, [])

    for single_file in file_clipboard:
        if os.path.exists(single_file):
            auth.create_data(project=project, path=single_file, user=request.user)

    clipboard[const.FILES_CLIPBOARD] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@read_access(type=Data)
def data_view(request, uid):
    "Show information specific to each data."

    data = Data.objects.get_all(uid=uid).first()
    project = data.project

    context = dict(data=data, project=project, activate='Selected Data')
    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@write_access(type=Data, fallback_view="data_view")
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
    context = dict(data=data, form=form)
    return render(request, 'data_edit.html', context)


@write_access(type=Project, fallback_view="data_list")
def data_upload(request, uid):
    "Data upload view routed through auth.create_data."

    owner = request.user
    project = Project.objects.get_all(uid=uid).first()
    form = forms.DataUploadForm(user=owner, project=project)

    if request.method == "POST":

        form = forms.DataUploadForm(data=request.POST, files=request.FILES, user=owner, project=project)

        if form.is_valid():
            data = form.save()
            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", request=request, kwargs={'uid': project.uid}))

    uploaded_files = Data.objects.filter(owner=owner, method=Data.UPLOAD)

    # The current size of the existing data
    current_size = uploaded_files.aggregate(Sum("size"))["size__sum"] or 0

    # Maximum data that may be uploaded.
    maximum_size = owner.profile.max_upload_size * 1024 * 1024

    context = dict(project=project, form=form, activate="Add Data", maximum_size=maximum_size,
                   current_size=current_size)

    counts = get_counts(project)

    context.update(counts)

    return render(request, 'data_upload.html', context)


@read_access(type=Analysis)
def recipe_view(request, uid):
    """
    Returns a recipe view based on its id.
    """
    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project
    context = dict(recipe=recipe, project=project, activate='Recipe View')

    # How many results for this recipe
    rcount = Job.objects.filter(analysis=recipe).count()
    counts = get_counts(project)

    try:
        # Fill in the script with json data.
        json_data = auth.fill_data_by_name(project=project, json_data=recipe.json_data)
        ctx = Context(json_data)
        script_template = Template(recipe.template)
        script = script_template.render(ctx)
    except Exception as exc:
        messages.error(request, f"Error rendering code: {exc}")
        script = recipe.template

    context.update(counts, rcount=rcount, script=script)

    return render(request, "recipe_view.html", context)


@read_access(type=Analysis)
def recipe_code_download(request, uid):
    """
    Download the raw recipe template as a file
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    try:
        # Fill in the script with json data.
        context = Context(recipe.json_data)
        script_template = Template(recipe.template)
        script = script_template.render(context)
    except Exception as exc:
        logger.error(exc)

    # Trigger file download with name of the recipe
    filename = "_".join(recipe.name.split()) + ".sh"

    response = HttpResponse(script, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename={filename}'

    return response


@read_access(type=Analysis)
def recipe_code_view(request, uid):
    """
    Returns an analysis code view based on its id.
    """
    user = request.user
    recipe = Analysis.objects.get_all(uid=uid).first()

    project = recipe.project

    try:
        # Fill in the script with json data.
        json_data = auth.fill_data_by_name(project=project, json_data=recipe.json_data)
        ctx = Context(json_data)
        script_template = Template(recipe.template)
        script = script_template.render(ctx)
    except Exception as exc:
        messages.error(request, f"Error rendering code: {exc}")
        script = recipe.template

    context = dict(recipe=recipe, project=project, activate='Recipe Code', script=script)

    rcount = Job.objects.filter(analysis=recipe).count()
    counts = get_counts(project)
    context.update(counts, rcount=rcount)

    return render(request, "recipe_code_view.html", context)


@read_access(type=Analysis)
@ratelimit(key='ip', rate='10/h', block=True, method=ratelimit.UNSAFE)
def recipe_run(request, uid):
    """
    View used to execute recipes and start a 'Queued' job.
    """

    analysis = Analysis.objects.get_all(uid=uid).first()

    project = analysis.project

    # Form submission.
    if request.method == "POST":

        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data,
                                     data=request.POST)

        # The form validation will authorize the job.
        if form.is_valid():

            # The desired name of for the results.
            name = form.cleaned_data.get("name")

            # Generates the JSON data from the bound form field.
            json_data = form.fill_json_data()

            # Create the job from the recipe and incoming json data.
            job = auth.create_job(analysis=analysis, user=request.user,
                                  json_data=json_data, name=name)

            # Spool the job right away if UWSGI exists.
            if tasks.HAS_UWSGI:
                # Update the job state.
                Job.objects.get_all(id=job.id).update(state=Job.SPOOLED)

                # Spool via UWSGI.
                tasks.execute_job.spool(job_id=job.id)

            return redirect(reverse("job_list", request=request, kwargs=dict(uid=project.uid)))
    else:
        initial = dict(name=f"Results for: {analysis.name}")
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, form=form, activate='Run Recipe')

    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


@read_access(type=Analysis)
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
        form = forms.EditCode(user=user, project=project, data=request.POST, recipe=analysis)

        if form.is_valid():
            # Preview action will let the form cascade through.
            commit = form.cleaned_data['action'] == 'SAVE'
            analysis = form.save(commit=commit)
            if commit:
                messages.info(request, "The recipe has been updated.")
                return redirect(reverse("recipe_view", request=request, kwargs=dict(uid=analysis.uid)))
    else:
        # This gets triggered on a GET request.
        initial = dict(template=analysis.template, json=analysis.json_text)
        form = forms.EditCode(user=user, project=project, initial=initial, recipe=analysis)

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


@write_access(type=Analysis, fallback_view='recipe_view')
def recipe_edit(request, uid):
    "Edit meta-data associated with a recipe."

    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project

    action_url = reverse('recipe_edit', request=request, kwargs=dict(uid=recipe.uid))
    form = forms.RecipeForm(instance=recipe, user=request.user)

    if request.method == "POST":
        form = forms.RecipeForm(data=request.POST, files=request.FILES, instance=recipe, user=request.user)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", request=request, kwargs=dict(uid=recipe.uid)))

    context = dict(analysis=recipe, project=project, form=form, action_url=action_url,
                   name=recipe.name)

    return render(request, 'recipe_edit.html', context)


@write_access(type=Job, fallback_view="job_view")
def job_edit(request, uid):
    "Edit meta-data associated with a job."

    job = Job.objects.get_all(uid=uid).first()
    project = job.project
    form = forms.JobEditForm(instance=job, user=request.user)

    if request.method == "POST":
        form = forms.JobEditForm(data=request.POST, files=request.FILES, instance=job, user=request.user)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", request=request, kwargs=dict(uid=job.uid)))

    context = dict(job=job, project=project, form=form)
    return render(request, 'job_edit.html', context)


@write_access(type=Analysis, fallback_view="recipe_view")
def recipe_delete(request, uid):
    recipe = Analysis.objects.get_all(uid=uid).first()

    auth.delete_object(obj=recipe, request=request)

    return redirect(reverse("recipe_list", kwargs=dict(uid=recipe.project.uid)))


@write_access(type=Job, fallback_view="job_view")
def job_delete(request, uid):
    job = Job.objects.get_all(uid=uid).first()

    running_job = job.state == Job.RUNNING and not job.deleted

    if running_job:
        messages.error(request, "Can not delete a running job. Wait until it finishes.")
        return redirect(job.url())

    auth.delete_object(obj=job, request=request)
    return redirect(reverse("job_list", kwargs=dict(uid=job.project.uid)))


@write_access(type=Data, fallback_view="data_view")
def data_delete(request, uid):
    data = Data.objects.get_all(uid=uid).first()

    auth.delete_object(obj=data, request=request)

    return redirect(reverse("data_list", kwargs=dict(uid=data.project.uid)))


@read_access(type=Job)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.get_all(uid=uid).first()
    project = job.project

    # The path is a GET parameter
    path = request.GET.get('path', "")

    context = dict(job=job, project=project, activate='View Result', path=path)

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


@read_access(type=Data)
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
