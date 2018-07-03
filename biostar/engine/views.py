import os
import logging

from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.db.models import Q

from sendfile import sendfile
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.utils import timezone
from django.urls import reverse
from django.utils.safestring import mark_safe
from biostar.accounts.models import Profile, User
from biostar.forum.models import Post

from . import tasks, auth, forms, util

from .diffs import color_diffs
from .decorators import object_access
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


@object_access(type=Project, access=Access.READ_ACCESS, url='project_view')
def clear_clipboard(request, uid, url="project_view", board=""):
    "Clear copied objects held in clipboard."

    if board:
        request.session[board] = None

    return redirect(reverse(url, kwargs=dict(uid=uid)))


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
            return redirect(reverse("project_users", kwargs=dict(uid=project.uid)))

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


def paste(project, post_request, board):
    "Used to paste data and job"

    form = forms.PasteForm(project=project, data=post_request.POST, request=post_request, board=board)
    if form.is_valid():
        form.save()
        post_request.session[board] = None
        return True, form
    else:
        return False, form


@object_access(type=Project, access=Access.READ_ACCESS)
def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """
    project = Project.objects.filter(uid=uid).first()
    if request.method == 'POST':
        success, form = paste(project=project, post_request=request, board="files_clipboard")
        if success:
            return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))
    else:
        form = forms.PasteForm(project=project, request=request, board='files_clipboard')

    return project_view(request=request, uid=uid, template_name="data_list.html",
                        active='data', extra_context=dict(form=form))


@object_access(type=Project, access=Access.READ_ACCESS)
def discussion_list(request, uid):

    project = Project.objects.filter(uid=uid).first()
    posts = project.post_set

    context = dict(posts=posts)
    return project_view(request=request, uid=uid, template_name="discussion_list.html",
                        active='discussion', extra_context=context)


@object_access(type=Post, access=Access.READ_ACCESS)
def discussion_view(request, uid):


    return


@object_access(type=Project, access=Access.READ_ACCESS)
def recipe_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    project = Project.objects.filter(uid=uid).first()
    if request.method == 'POST':
        success, form = paste(project=project, post_request=request, board="recipe_clipboard")
        if success:
            return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))
    else:
        form = forms.PasteForm(project=project, request=request, board='recipe_clipboard')

    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes',
                        more_info=True, extra_context=dict(form=form))


def job_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


def get_counts(project):
    data_count = Data.objects.filter(project=project).count()
    recipe_count = Analysis.objects.filter(project=project).count()
    result_count = Job.objects.filter(project=project).count()
    discussion_count = Post.objects.exclude(status=Post.DELETED).filter(project=project).count()
    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count,
        discussion_count=discussion_count
    )


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, uid, template_name="recipe_list.html", active='recipes', more_info=None,
                 extra_context={}):

    project = Project.objects.filter(uid=uid).first()
    # Show counts for the project.
    counts = get_counts(project)

    # Select all the data in the project.
    data_list = Data.objects.filter(project=project).order_by("-sticky", "-date").all()
    recipe_list = Analysis.objects.filter(project=project).order_by("-sticky", "-date").all()
    job_list = Job.objects.filter(project=project).order_by("-sticky", "-date").all()

    # Filter job results by analysis
    filter = request.GET.get('filter', '')
    if filter:
        filter = Analysis.objects.filter(uid=filter).first()
        job_list = job_list.filter(analysis=filter)

    # This is not quite right to be here.
    if more_info:
        more_info = project.html

    context = dict(project=project, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, filter=filter, more_info=more_info)
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
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

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
            name = form.cleaned_data["name"]
            text = form.cleaned_data["text"]
            summary = form.cleaned_data["summary"]
            stream = form.cleaned_data["image"]
            sticky = form.cleaned_data["sticky"]
            privacy = form.cleaned_data["privacy"]
            owner = request.user
            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream, sticky=sticky, privacy=privacy)
            project.save()
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

    context = dict(form=form)
    return render(request, "project_create.html", context=context)


@object_access(type=Data, access=Access.READ_ACCESS)
def data_view(request, uid):
    "Show information specific to each data."

    data = Data.objects.get_all(uid=uid).first()
    project = data.project
    path = request.GET.get('path', '')

    if request.method == "POST":
        form = forms.FileCopyForm(request, data.uid, data.get_data_dir(), request.POST)
        if form.is_valid():
            # Copies data to clipboard
            ndata = form.save()
            msg = mark_safe(f"Copied <b>{ndata}</b> files from <b>{project.name}</b> to the Clipboard.")
            messages.success(request, msg)
            return redirect(reverse("data_view", kwargs=dict(uid=data.uid)))
    else:
        form = forms.FileCopyForm(request, data.uid, data.get_data_dir())

    root = data.get_data_dir()
    abspath = join(root, path)
    try:
        files = util.scan_files(abspath=abspath, relpath=path, root=root)
    except Exception as exc:
        messages.error(request, f"{exc}")
        files = []

    context = dict(data=data, project=project, activate='Selected Data', files=files, path=path, form=form)
    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.OWNER_ACCESS, url='data_view')
def data_edit(request, uid):
    """
    Edit meta-data associated with Data.
    """

    data = Data.objects.filter(uid=uid).first()
    form = forms.DataEditForm(instance=data, initial=dict(type=data.type), user=request.user)

    if request.method == "POST":
        form = forms.DataEditForm(data=request.POST, instance=data, user=request.user, files=request.FILES)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(uid=data.uid)))
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
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))

    context = dict(project=project, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER)
def recipe_view(request, uid):
    """
    Returns an analysis view based on its id.
    """
    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project

    if request.method == "POST":
        form = forms.RecipeCopyForm(data=request.POST, recipe=recipe, request=request)
        if form.is_valid():
            form.save()
            return redirect(reverse("recipe_view", kwargs=dict(uid=recipe.uid)))
    else:
        form = forms.RecipeCopyForm(recipe=recipe, request=request)

    context = dict(recipe=recipe, project=project, activate='Selected Recipe', form=form)
    counts = get_counts(project)
    context.update(counts)
    return render(request, "recipe_view.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view')
def recipe_run(request, uid):
    """
    View used to execute recipes and start a 'Queued' job.
    """

    analysis = Analysis.objects.filter(uid=uid).first()
    project = analysis.project

    if request.method == "POST":
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, data=request.POST)

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

            return redirect(reverse("job_list", kwargs=dict(uid=project.uid)))
    else:
        initial = dict(name=f"Results for: {analysis.name}")
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, form=form, activate='Run Recipe')
    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER, url='recipe_view')
def recipe_code(request, uid):
    """
    Displays and allows edit on a recipe code.

    Since we allow a preview for un-authenticated users thus the view
    is more complicated than a typical DJANGO form handler.
    """
    user = request.user

    # There has to be a recipe to work with.
    analysis = Analysis.objects.filter(uid=uid).first()
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
                return redirect(reverse("recipe_view", kwargs=dict(uid=analysis.uid)))
    else:
        # This gets triggered on a GET request.
        initial = dict(template=analysis.template, json=analysis.json_text)
        form = forms.EditCode(user=user, project=project, initial=initial)

    # Bind the JSON to the form.
    recipe = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=dict(name=name))

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
    return render(request, 'recipe_code.html', context)


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

            return redirect(reverse('recipe_list', kwargs=dict(uid=project.uid)))
    # The url to submit to.
    action_url = reverse('recipe_create', kwargs=dict(uid=project.uid))
    context = dict(project=project, form=form, action_url=action_url, name="New Recipe")

    return render(request, 'recipe_edit.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MANAGER,url='recipe_view')
def recipe_diff(request, uid):
    """
    View used to show diff in template and authorize it.
    Restricted to moderators and staff members.
    """
    recipe = Analysis.objects.filter(uid=uid).first()
    differ = auth.template_changed(template=recipe.last_valid, analysis=recipe)
    differ = color_diffs(differ)

    form = forms.RecipeDiff(recipe=recipe, request=request, user=request.user)

    if request.method == "POST":
        form = forms.RecipeDiff(recipe=recipe, user=request.user, data=request.POST,
                          request=request)
        if form.is_valid():
            form.save()
            return redirect(reverse('recipe_view', kwargs=dict(uid=recipe.uid)))

    context = dict(activate="Recent Template Change",  project=recipe.project, recipe=recipe,
                   diff=mark_safe(''.join(differ)), form=form)

    counts = get_counts(recipe.project)
    context.update(counts)

    return render(request, "recipe_diff.html", context=context)



@object_access(type=Analysis, access=Access.OWNER_ACCESS, url='recipe_view')
def recipe_edit(request, uid):
    "Edit meta-data associated with a recipe."

    recipe = Analysis.objects.get_all(uid=uid).first()
    project = recipe.project
    action_url = reverse('recipe_edit', kwargs=dict(uid=recipe.uid))
    form = forms.RecipeForm(instance=recipe)

    if request.method == "POST":
        form = forms.RecipeForm(data=request.POST, files=request.FILES, instance=recipe)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", kwargs=dict(uid=recipe.uid)))

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
            return redirect(reverse("job_view", kwargs=dict(uid=job.uid)))

    context = dict(job=job, project=project, form=form)
    return render(request, 'job_edit.html', context)


def object_state_toggle(request, uid, obj_type):
    """
    Toggle instance.deleted if user has owner access to instance.
    """

    # Map obj_type to an object and respective url
    obj_map = dict(job=(Job, 'job_list'),
                   data=(Data, 'data_list'),
                   recipe=(Analysis, "recipe_list"))
    if not obj_map.get(obj_type):
        messages.error(request, "Can not toggle state.")
        return redirect(reverse('project_list'))

    # Make query and build urls
    obj, view_name = obj_map[obj_type][0], obj_map[obj_type][1]
    instance = obj.objects.get_all(uid=uid).first()

    # Ensure running jobs do not get deleted.
    if isinstance(instance, Job) and instance.state == Job.RUNNING and not instance.deleted:
        messages.error(request, "Can not delete a running job. Wait until it finishes.")
        return redirect(instance.url())

    # Make sure user has owner access to instance before toggling
    has_access = auth.check_obj_access(instance=instance, user=request.user, request=request,
                                       access=Access.OWNER_ACCESS)
    msg = f"Deleted <b>{instance.name}</b>. View in Recycle Bin."
    name_repeat = auth.check_data_name(name=instance.name, data=instance, bool=True)

    if has_access:
        if name_repeat and instance.deleted and isinstance(instance, Data):
            messages.error(request, "Name already exists. Edit the data then restore.")
            return redirect(instance.url())

        # Toggle delete state
        instance.deleted = not instance.deleted
        instance.save()
        msg = msg if instance.deleted else f"Restored <b>{instance.name}</b>."
        messages.success(request, mark_safe(msg))

    return redirect(reverse(view_name, kwargs=dict(uid=instance.project.uid)))



@object_access(type=Job, access=Access.READ_ACCESS)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.get_all(uid=uid).first()
    project = job.project

    # The path is a GET parameter
    path = request.GET.get('path', "")

    if request.method == "POST":
        form = forms.FileCopyForm(request, job.uid, job.path, request.POST)
        if form.is_valid():
            # Copies data to clipboard
            ndata = form.save()
            msg = mark_safe(f"Copied <b>{ndata}</b> files from <b>{job.name}</b> to the Clipboard.")
            messages.success(request, msg)
            return redirect(reverse("job_list", kwargs=dict(uid=project.uid)))
    else:
        form = forms.FileCopyForm(request, job.uid, job.path)

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

    context = dict(job=job, project=project, files=files, path=path, activate='Selected Result', form=form)

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def file_serve(request, uid, path, obj):
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
    return file_serve(request=request, path=path, uid=uid, obj=obj)


@object_access(type=Job, access=Access.READ_ACCESS, url='job_view')
def job_serve(request, uid, path):
    """
    Serves files from a job directory.
    """
    obj = Job.objects.get_all(uid=uid).first()
    return file_serve(request=request, path=path, uid=uid, obj=obj)
