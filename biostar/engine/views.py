import glob
import logging
import mistune

from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.db.models import Q

from sendfile import sendfile
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse
from django.utils.safestring import mark_safe
from biostar.accounts.models import Profile


from . import tasks, util
from .decorators import object_access
from .forms import *
from .models import (Project, Data, Analysis, Job, Access)

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
__DOCS_DIR = join(__CURRENT_DIR, "docs")

logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


def docs(request, name):
    patt = join(__DOCS_DIR, name) + ".*"
    files = glob.glob(patt)
    if not files:
        msg = f"Cannot be find the requested page: {name} "
        messages.error(request, msg)
        return redirect("index")
    if len(files) > 1:
        msg = f"Multiple files match: {{name}}"
        messages.warning(request, msg)
    target = files[0]
    content = open(target).read()

    # Render markdown into HTML.
    if target.endswith(".md"):
        content = mistune.markdown(content)

    title = name.replace("-", " ").replace("_", " ").title()
    context = dict(content=content, title=title)
    return render(request, 'info/doc_base.html', context=context)


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


def toggle_notifications(request):
    "Allows user to toggle email notification to projects they have access to "

    form = ToggleNotifications(user=request.user)

    if request.method == "POST":
        form = ToggleNotifications(data=request.POST, user=request.user)
        if form.is_valid():
            form.save()
            return redirect(reverse("profile"))


    user = User.objects.filter(id=id).first()
    context = dict(user=user, form=form)

    return render(request, 'accounts/profile.html', context)



def recycle_bin(request):
    "Recycle bin view for a user"

    if request.user.is_anonymous:
        messages.error(request, "You must be logged in to view recycle bin.")
        return redirect("/")

    # Only searches projects you have access.
    all_projects = auth.get_project_list(user=request.user)

    del_data = Data.objects.filter(deleted=True, project__in=all_projects,
                                   owner=request.user).order_by("date")

    del_recipes = Analysis.objects.filter(deleted=True, project__in=all_projects,
                                          owner=request.user).order_by("date")

    del_jobs = Job.objects.filter(deleted=True, project__in=all_projects,
                                  owner=request.user).order_by("date")

    context = dict(jobs=del_jobs, data=del_data, recipes=del_recipes)

    return render(request, 'recycle_bin.html', context=context)


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
    form = ChangeUserAccess()

    if request.method == "POST":
        form = ChangeUserAccess(data=request.POST)

        # User needs to be authenticated and have admin access to make any changes.
        if form.is_valid() and request.user.is_authenticated:
            user, access = form.save()
            msg = f"Changed <b>{user.first_name}</b>'s access to {label(access.get_access_display())}"
            messages.success(request, mark_safe(msg))
            return redirect(reverse("project_users", kwargs=dict(uid=project.uid)))

    # Users that have been searched for.
    targets = User.objects.filter(Q(email__contains=q) | Q(first_name__contains=q)) if q else []
    current = access_forms(users=user_list, project=project, exclude=[request.user])
    results = access_forms(users=targets, project=project, exclude=[request.user])
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


def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """

    return project_view(request=request, uid=uid, template_name="data_list.html", active='data')


def recipe_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes', more_info=True)


def job_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


@object_access(type=Project, access=Access.READ_ACCESS, url='data_list')
def data_navigate(request, uid):
    """
    Renders project data as if they were files.
    """

    project = Project.objects.filter(uid=uid).first()
    # Same ordering as data_list
    all_data = project.data_set.order_by("sticky", "-date").all()

    if request.method == "POST":
        form = DataCopyForm(data=request.POST, request=request)
        if form.is_valid():
            # Copies data to clipboard
            ndata = form.save()
            msg = mark_safe(f"Copied <b>{ndata}</b> data from <b>{project.name}</b> to the Clipboard.")
            messages.success(request, msg)
            return redirect(reverse("data_navigate", kwargs=dict(uid=project.uid)))
    else:
        form = DataCopyForm(request=request)

    context = dict(activate='File Explorer', all_data=all_data, project=project, form=form)
    counts = get_counts(project)
    context.update(counts)
    return render(request, "data_navigator.html", context)


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_browse(request, uid, path=''):
    """
    Browse the directory that corresponds to job results.
    """
    job = Job.objects.filter(uid=uid).first()
    project = job.project

    form = FileCopyForm(root_dir=job.path, request=request, uid=job.uid)
    if request.method == "POST":
        form = FileCopyForm(data=request.POST, uid=job.uid, request=request, root_dir=job.path)
        if form.is_valid():
            # Copies files to clipboard
            nfiles = form.save()
            msg = mark_safe(f"Copied <b>{nfiles}</b> file(s) to Clipboard from <b>{job.name}</b>")
            messages.success(request, msg)
            return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))

    context = dict(activate='Selected Job', job=job, project=project, form=form)
    counts = get_counts(project)
    context.update(counts)

    return file_browse(request=request, instance=job, path=path,
                       template_name="job_browser.html", extra_context=context)


@object_access(type=Data, access=Access.READ_ACCESS, url='data_view')
def data_browse(request, uid, path=''):
    """
    Browse the directory that corresponds to data.
    """

    data = Data.objects.filter(uid=uid).first()
    project = data.project

    form = FileCopyForm(root_dir=data.get_data_dir(), request=request, uid=data.uid)
    if request.method == "POST":
        form = FileCopyForm(data=request.POST, uid=data.uid, request=request, root_dir=data.get_data_dir())
        if form.is_valid():
            # Copies files to clipboard
            nfiles = form.save()
            msg = mark_safe(f"Copied <b>{nfiles}</b> file(s) to Clipboard from <b>{data.name}</b>")
            messages.success(request, msg)
            return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))

    back_uid = None if path else project.uid
    context = dict(activate='File Explorer', data=data, project=project, project_uid=back_uid, form=form)

    counts = get_counts(project)
    context.update(counts)

    return file_browse(request=request, instance=data, path=path,
                       template_name="data_browser.html", extra_context=context)


def file_browse(request, instance, template_name, path='', extra_context={}):
    """
    Browse the filesystem that corresponds to a root folder.
    """

    # Instance is expected to be a Job or Data object.
    if isinstance(instance, Job):
        root, exclude = instance.path, ''
    else:
        # Exclude toc from file_list
        exclude = os.path.basename(instance.get_path())
        root = instance.get_data_dir()

    # Get the root directory.
    target_path = join(root, path)

    if target_path.startswith(root) and os.path.exists(target_path):
        # Pathlike objects with attributes such as name, is_file
        file_list = list(filter(lambda p: p.name != exclude, os.scandir(target_path)))
        # Sort the file list. Directories first, then by name.
        file_list = sorted(file_list, key=lambda p: (p.is_file(), p.name))
    else:
        # Attempting to access a file outside of the root directory
        messages.error(request, "Invalid path.")
        file_list = []

    context = dict(file_list=file_list, instance=instance, path=path)
    context.update(extra_context)

    return render(request, template_name, context)


def get_counts(project):
    data_count = Data.objects.filter(deleted=False, project=project).count()
    recipe_count = Analysis.objects.filter(deleted=False, project=project).count()
    result_count = Job.objects.filter(deleted=False, project=project).count()
    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count
    )


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, uid, template_name="recipe_list.html", active='recipes', more_info=None):
    project = Project.objects.filter(uid=uid).first()
    # Show counts for the project.
    counts = get_counts(project)

    # Select all the data in the project.
    data_list = Data.objects.filter(deleted=False, project=project).order_by("sticky", "-date").all()
    recipe_list = Analysis.objects.filter(deleted=False, project=project).order_by("-date").all()
    job_list = Job.objects.filter(deleted=False, project=project).order_by("-date").all()

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

    return render(request, template_name, context)


@object_access(type=Project, access=Access.OWNER_ACCESS, url='data_list')
def project_edit(request, uid):
    "Edit meta-data associated with a project."

    project = Project.objects.filter(uid=uid).first()
    form = ProjectForm(instance=project)

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES, instance=project)
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
    form = ProjectForm(initial=initial)

    if request.method == "POST":
        # create new projects here ( just populates metadata ).
        form = ProjectForm(request.POST, request.FILES)
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

    data = Data.objects.filter(uid=uid).first()

    project = data.project
    context = dict(data=data, project=project, activate='Selected Data')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.READ_ACCESS, url='data_view')
def data_copy(request, uid):
    "Store Data object in data clipboard "

    data = Data.objects.filter(uid=uid).first()

    project = data.project

    auth.load_data_clipboard(uid=uid, request=request)

    messages.success(request, mark_safe(f"Copied data <b>{data.name}</b> to the clipboard."))

    return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.WRITE_ACCESS, url='data_list')
def data_paste(request, uid):
    """
    Paste data in clipboard to project.
    """

    project = Project.objects.filter(uid=uid).first()

    # A list of data objects in the data_clipboard
    data_set = auth.dump_data_clipboard(request, reset=True)

    for data in data_set:
        # Create data in project by linking files ( excluding toc file ).
        new_data = auth.create_data(project=project, name=f"Copy of {data.name}", text=data.text,
                                    path=data.get_data_dir(), summary=data.summary,
                                    type=data.type, skip=data.get_path(), user=request.user)
        new_data.state = data.state
        new_data.save()

    if data_set:
        msg = mark_safe(f"Pasted <b>{len(data_set)}</b> data to <b>{project.name}</b>.")
        messages.success(request, msg)

    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.WRITE_ACCESS, url='data_list')
def files_paste(request, uid):
    """
    View used to paste result files copied from a job or data.
    """

    files = request.session.get("files_clipboard", None)
    project = Project.objects.filter(uid=uid).first()
    url = reverse("data_list", kwargs=dict(uid=project.uid))

    root_path = auth.validate_files_clipboard(clipboard=files, request=request)
    if not root_path:
        return redirect(url)

    # Some files in clipboard might be outside job path.
    files = [f for f in files if f.startswith(root_path)]
    # Add data to project
    for file in files:
        auth.create_data(project=project, path=file, user=request.user)

    request.session["files_clipboard"] = None
    msg = mark_safe(f"Pasted <b>{len(files)}</b> file(s) to project <b>{project.name}</b>.")
    messages.success(request, msg)
    return redirect(url)


@object_access(type=Data, access=Access.OWNER_ACCESS, url='data_view')
def data_edit(request, uid):
    """
    Edit meta-data associated with Data.
    """

    data = Data.objects.filter(uid=uid).first()
    form = DataEditForm(instance=data, initial=dict(type=data.type))

    if request.method == "POST":
        form = DataEditForm(data=request.POST, instance=data)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(uid=data.uid)))

    context = dict(data=data, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.WRITE_ACCESS, url='data_list')
def data_upload(request, uid):
    "Data upload view routed through auth.create_data."

    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    form = DataUploadForm(user=owner)
    if request.method == "POST":
        form = DataUploadForm(data=request.POST, files=request.FILES, user=owner)

        if form.is_valid():
            text = form.cleaned_data["text"]
            stream = form.cleaned_data["file"]
            summary = form.cleaned_data["summary"]
            type = form.cleaned_data["type"]
            name = stream.name
            data = auth.create_data(stream=stream, name=name,
                                    text=text, user=owner, project=project, summary=summary,
                                    type=type)

            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))

    context = dict(project=project, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
def recipe_view(request, uid):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(uid=uid).first()
    project = analysis.project
    context = dict(analysis=analysis, project=project, activate='Selected Recipe')

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
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, data=request.POST)

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
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, form=form, activate='Run Recipe')
    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view')
def recipe_copy(request, uid):
    """
    Store Analysis object in request.sessions['recipe_clipboard']
    """

    recipe = Analysis.objects.filter(uid=uid).first()
    project = recipe.project
    request.session["recipe_clipboard"] = recipe.uid

    msg = f"Copied recipe <b>{recipe.name}</b> to Clipboard."
    msg = mark_safe(msg)
    messages.success(request, msg)
    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.WRITE_ACCESS, url='project_view')
def recipe_paste(request, uid):
    """
    Paste recipe stored in the clipboard into project.
    """

    recipe_uid = request.session.get("recipe_clipboard")
    recipe = Analysis.objects.filter(uid=recipe_uid).first()
    project = Project.objects.filter(uid=uid).first()
    url = reverse("recipe_list", kwargs=dict(uid=project.uid))

    # Make sure user has read access to recipe in clipboard
    has_access = auth.check_obj_access(instance=recipe, user=request.user, request=request,
                                       access=Access.READ_ACCESS)
    if not has_access:
        request.session["recipe_clipboard"] = None
        return redirect(url)

    attrs = auth.get_analysis_attr(recipe, project=project)
    attrs.update(stream=recipe.image, name=f"Copy of {recipe.name}", security=recipe.security,
                 user=request.user)
    new_recipe = auth.create_analysis(**attrs)
    # Ensure the diff gets inherited.
    new_recipe.last_valid = recipe.last_valid
    new_recipe.save()

    msg = f"Pasted recipe <b>{recipe.name}</b> to project <b>{project.name}</b>."
    msg = mark_safe(msg)
    messages.success(request, msg)
    request.session["recipe_clipboard"] = None

    return redirect(url)


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MODERATOR, url='recipe_view')
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
        form = EditCode(user=user, project=project, data=request.POST)

        if form.is_valid():

            # Templates.
            template = form.cleaned_data['template']

            # Preview action will let the form cascade through.
            save = form.cleaned_data['action'] == 'SAVE'

            analysis.json_text = form.cleaned_data['json']

            # Changes to template will require a review ( only when saving ).
            if auth.template_changed(analysis=analysis, template=template) and save:
                analysis.security = Analysis.UNDER_REVIEW

            # Moderators and staff members will automatically get authorized.
            if user.is_staff:# or user.profile.is_moderator:
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
        form = EditCode(user=user, project=project, initial=initial)

    # Bind the JSON to the form.
    recipe = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=dict(name=name))

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
    form = RecipeForm(initial=dict(name="New Recipe"))

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES)

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


@object_access(type=Analysis, access=Access.READ_ACCESS, role=Profile.MODERATOR, url='recipe_view')
def recipe_diff(request, uid):
    """
    View used to show diff in template and authorize it.
    Restricted to moderators and staff members.
    """
    recipe = Analysis.objects.filter(uid=uid).first()

    differ = auth.template_changed(template=recipe.last_valid, analysis=recipe)
    differ = auth.color_diffs(differ)
    context = dict(activate="Recent Template Change",  project=recipe.project, recipe=recipe,
                   diff=mark_safe(''.join(differ)))
    counts = get_counts(recipe.project)
    context.update(counts)

    return render(request, "recipe_diff.html", context=context)


# Ensure only moderators access when role=moderator and access=Access.NO_ACCESS
@object_access(type=Analysis, role=Profile.MODERATOR, access=Access.NO_ACCESS, url='recipe_diff')
def recipe_approve(request, uid):
    "Approve changes made in recipe. Only moderators are allowed this action."

    recipe = Analysis.objects.filter(uid=uid).first()

    recipe.last_valid = recipe.template
    recipe.security = Analysis.AUTHORIZED
    recipe.save()

    messages.success(request, "Recipe changes have been approved.")
    return redirect(reverse('recipe_diff', kwargs=dict(uid=recipe.uid)))


@object_access(type=Analysis, access=Access.WRITE_ACCESS, role=Profile.MODERATOR, url='recipe_diff')
def recipe_revert(request, uid):
    """
        Allowed to moderators and users with write access.
        Revert changes made to recipes back to original.
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    recipe.template = recipe.last_valid
    recipe.security = Analysis.AUTHORIZED
    recipe.save()

    messages.success(request, "Recipe has been reverted to original.")
    return redirect(reverse('recipe_diff', kwargs=dict(uid=recipe.uid)))


@object_access(type=Analysis, access=Access.OWNER_ACCESS, url='recipe_view')
def recipe_edit(request, uid):
    "Edit meta-data associated with a recipe."

    recipe = Analysis.objects.filter(uid=uid).first()
    project = recipe.project
    action_url = reverse('recipe_edit', kwargs=dict(uid=recipe.uid))
    form = RecipeForm(instance=recipe)

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES, instance=recipe)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", kwargs=dict(uid=recipe.uid)))

    context = dict(analysis=recipe, project=project, form=form, action_url=action_url,
                   name=recipe.name)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Job, access=Access.OWNER_ACCESS, url="job_view")
def job_edit(request, uid):
    "Edit meta-data associated with a job."

    job = Job.objects.filter(uid=uid).first()
    project = job.project
    form = JobEditForm(instance=job)

    if request.method == "POST":
        form = JobEditForm(data=request.POST, files=request.FILES, instance=job)
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
    instance = obj.objects.filter(uid=uid).first()
    url = reverse(view_name, kwargs=dict(uid=instance.project.uid))
    bin_url = reverse('recycle_bin')
    # Make sure user has owner access to instance before toggling
    has_access = auth.check_obj_access(instance=instance, user=request.user, request=request,
                                       access=Access.OWNER_ACCESS)

    msg = f"Deleted <b>{instance.name}</b>. View in <a href={bin_url}>Recycle Bin</a>."
    if has_access:
        # Toggle delete state
        instance.deleted = not instance.deleted
        instance.save()
        msg = msg if instance.deleted else f"Restored <b>{instance.name}</b>."
        messages.success(request, mark_safe(msg))

    return redirect(url)


@object_access(type=Job, access=Access.READ_ACCESS)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.filter(uid=uid).first()
    project = job.project

    context = dict(job=job, project=project, activate='Selected Result')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def file_serve(request, uid, path, klass):
    """
    Authenticates access through decorator before serving file.
    """

    # Get the object that corresponds to the entry.
    obj = klass.objects.filter(uid=uid).first()

    # Different file layout for jobs and data.
    if isinstance(obj, Data):
        url = 'data_entry'
        root = obj.get_data_dir()
    else:
        # Job type.
        root = obj.path
        url = 'job_entry'

    # This will turn into an absolute path.
    file_path = join(root, path)

    # This ensure only files in the object root can be accessed.
    if not file_path.startswith(root):
        messages.error(request, "Invalid path.")
        return redirect(reverse(url, kwargs=dict(uid=uid)))

    # Check that the file exists.
    if not os.path.isfile(file_path):
        messages.error(request, f"File not found: {path}")
        return redirect(reverse(url, kwargs=dict(uid=uid)))

    # The response will be the file content.
    mimetype = auth.guess_mimetype(fname=path)

    # Get the filesize in Mb
    size = os.path.getsize(file_path)/1024/1024

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
    return file_serve(request=request, path=path, uid=uid, klass=Data)


@object_access(type=Job, access=Access.READ_ACCESS, url='job_entry')
def job_serve(request, uid, path):
    """
    Serves files from a job directory.
    """
    return file_serve(request=request, path=path, uid=uid, klass=Job)
