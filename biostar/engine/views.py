import glob
import logging

from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.db.models import Q
import mistune

from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse

from . import tasks
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


def recycle_bin(request):
    "Recycle bin view for a user"

    if request.user.is_anonymous:
        messages.error(request, "You must be logged in to view recycle bin.")
        return redirect("/")

    all_projects = auth.get_project_list(user=request.user)

    # Only jobs enabled currently
    del_jobs = Job.objects.filter(state=Job.DELETED, project__in=all_projects,
                                  owner=request.user).order_by("date")

    context = dict(jobs=del_jobs, projects=[],recipe=[], data=[])

    return render(request, 'recycle_bin.html', context=context)


@object_access(type=Project, access=Access.READ_ACCESS, url='project_view')
def clear_clipboard(request, uid, redir="project_view", board=""):
    "Clear copy object held in clipboard"
    clear = [""] if board == "files_clipboard" else None

    if board:
        request.session[board] = clear

    return redirect(reverse(redir, kwargs=dict(uid=uid)))


def get_access(request, project):

    user = request.user if request.user.is_authenticated else None
    user_access = Access.objects.filter(project=project, user=user).first()
    # Current users access
    user_access = user_access or Access(access=Access.NO_ACCESS)

    # Users already with access to current project
    user_list = [a.user for a in project.access_set.all() if a.access > Access.NO_ACCESS]

    return user_access, user_list


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='data_list')
def project_users(request, uid):
    """
    Manage project users
    """
    project = Project.objects.filter(uid=uid).first()
    user_access, user_list = get_access(request, project)
    label = lambda x: f"<span class='ui green tiny label'>{x}</span>"

    # Search query
    q = request.GET.get("q", "")
    form = ChangeUserAccess()

    if request.method == "POST":
        form = ChangeUserAccess(data=request.POST)

        # User needs to be authenticated and have admin access to make any changes.
        if form.is_valid() and request.user.is_authenticated:
            user, access = form.change_access()
            msg = f"Changed <b>{user.first_name}</b>'s access to {label(access.get_access_display())}"
            messages.success(request, mark_safe(msg))
            return redirect(reverse("project_users", kwargs=dict(uid=project.uid)))

    # Users that have been searched for.
    targets = User.objects.filter(Q(email__contains=q) | Q(first_name__contains=q)) if q else []
    current = access_forms(users=user_list, project=project)
    results = access_forms(users=targets, project=project)
    context = dict(current=current, project=project, results=results, form=form, activate='selection',
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
    Returns the list of data for a project id.
    """

    return project_view(request=request, uid=uid, template_name="data_list.html", active='data')


def recipe_list(request, uid):
    """
    Returns the list of recipes for a project id.
    """
    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes', more_info=True)


def job_list(request, uid):
    """
    Returns the list of recipes for a project id.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


def files_list(request, instance, template_name, path='', extra_context={}):
    "File navigator used for  data and jobs"

    # Instance is expected to be a Job or Data object.
    if isinstance(instance, Job):
        root, exclude = instance.path, ''
    else:
        # Exclude toc from file_list
        exclude = os.path.basename(instance.get_path())
        root = instance.get_data_dir()

    target_path = join(root, path)

    if not target_path.startswith(root) or (not os.path.exists(target_path)):
        # Attempting to access a file outside of the root directory
        messages.error(request, "Path not in directory.")
        file_list = []
    else:
        # These are pathlike objects with attributes such as name, is_file
        file_list = list(filter(lambda p: p.name != exclude, os.scandir(target_path)))
        # Sort by properties
        file_list = sorted(file_list, key=lambda p: (p.is_file(), p.name))

    context = dict(file_list=file_list, instance=instance, path=path)
    context.update(extra_context)

    return render(request, template_name, context)


def get_counts(project):
    data_count = Data.objects.filter(project=project).count()
    recipe_count = Analysis.objects.filter(project=project).count()
    result_count = Job.objects.filter(~Q(state=Job.DELETED),project=project).count()
    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count
    )


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, uid, template_name="recipe_list.html", active='recipes', more_info=None):

    project = Project.objects.filter(uid=uid).first()
    # Show counts for the project.
    counts = get_counts(project)

    # Select all the data in the project.
    data_list = Data.objects.filter(project=project).order_by("sticky", "-date").all()
    recipe_list = Analysis.objects.filter(project=project).order_by("-date").all()
    job_list = Job.objects.filter(~Q(state=Job.DELETED), project=project).order_by("-date").all()

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


@object_access(type=Project, access=Access.EDIT_ACCESS, url='project_view', owner_only=True)
def project_edit(request, uid):

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
    data = Data.objects.filter(uid=uid).first()

    project = data.project
    context = dict(data=data, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.READ_ACCESS, url='data_view')
def data_copy(request, uid):
    "Store Data object in request.sessions['data_clipboard'] "

    data = Data.objects.filter(uid=uid).first()
    project = data.project
    request.session["data_clipboard"] = data.uid
    messages.success(request, mark_safe(f"Copied <b>{data.name}</b> to Clipboard."))

    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def data_paste(request, uid):
    "Paste data stored in the clipboard into project"

    data_uid = request.session.get("data_clipboard")
    data = Data.objects.filter(uid=data_uid).first()
    project = Project.objects.filter(uid=uid).first()

    # Make sure user has read access to data in clipboard
    has_access = auth.check_obj_access(instance=data, user=request.user, request=request,
                                       access=Access.READ_ACCESS)
    if not has_access:
        request.session["data_clipboard"] = None
        return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))

    # Create data object in project by linking files ( excluding toc file ).
    auth.create_data(project=project, name=f"Copy of {data.name}", text=data.text,
                     path=data.get_data_dir(),summary=data.summary,
                     type=data.type, skip=data.get_path())

    # Clear clipboard
    request.session["data_clipboard"] = None
    messages.success(request, mark_safe(f"Pasted <b>{data.name}</b> to <b>{project.name}</b>."))
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def files_paste(request, uid):
    "Paste files copied from a job to a project"

    files = request.session.get("files_clipboard", [""])
    project = Project.objects.filter(uid=uid).first()
    url = reverse("data_list", kwargs=dict(uid=project.uid))

    # Last item in clipboard is job.uid that the files in clipboard belong to.
    job_uid = files.pop(-1)
    job = Job.objects.filter(uid=job_uid).first()
    if not job:
        messages.error(request, "Files do not belong to any result")
        return redirect(url)

    # Make sure user has read access to job in clipboard.
    has_access = auth.check_obj_access(instance=job, user=request.user, request=request,
                                       access=Access.READ_ACCESS)
    if not has_access:
        request.session["files_clipboard"] = [""]
        return redirect(url)

    # Some files in clipboard might be outside job path.
    files = [f for f in files if f.startswith(job.path) ]
    # Add data to project
    for file in files:
        auth.create_data(project=project, path=file)

    request.session["files_clipboard"] = [""]
    msg = mark_safe(f"Pasted <b>{len(files)}</b> file(s) to project <b>{project.name}</b>.")
    messages.success(request, msg)
    return redirect(url)


@object_access(type=Data, access=Access.EDIT_ACCESS, url='data_view', owner_only=True)
def data_edit(request, uid):
    data = Data.objects.filter(uid=uid).first()
    form = DataEditForm(instance=data, initial=dict(type=data.type))

    if request.method == "POST":
        form = DataEditForm(data=request.POST, instance=data)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(uid=data.uid)))

    context = dict(data=data, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.READ_ACCESS, url='data_list')
def data_nav(request, uid):
    "Return special dir view of data list"

    project = Project.objects.filter(uid=uid).first()
    # Same ordering as data_list
    all_data = project.data_set.order_by("sticky", "-date").all()
    context = dict(activate='selection', all_data=all_data, project=project)
    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_nav.html", context)


@object_access(type=Project, access=Access.UPLOAD_ACCESS, url='data_list')
def data_upload(request, uid):
    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    form = DataUploadForm()
    if request.method == "POST":
        form = DataUploadForm(data=request.POST, files=request.FILES)

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


@object_access(type=Data, access=Access.READ_ACCESS, url='data_view')
def data_files_list(request, uid, path=''):
    data = Data.objects.filter(uid=uid).first()
    project = data.project

    back_uid = None if path else project.uid
    context = dict(activate='selection', data=data, project=project, project_uid=back_uid)

    counts = get_counts(project)
    context.update(counts)

    return files_list(request=request, instance=data, path=path,
                      template_name="data_files_list.html", extra_context=context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
def recipe_view(request, uid):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(uid=uid).first()
    project = analysis.project
    context = dict(analysis=analysis, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)
    return render(request, "recipe_view.html", context)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
def recipe_run(request, uid):
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

    context = dict(project=project, analysis=analysis, form=form, activate='selection')
    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view')
def recipe_copy(request, uid):
    "Store Analysis object in request.sessions['recipe_clipboard'] "

    recipe = Analysis.objects.filter(uid=uid).first()
    project = recipe.project
    request.session["recipe_clipboard"] = recipe.uid

    msg = f"Copied recipe <b>{recipe.name}</b> to Clipboard."
    msg = mark_safe(msg)
    messages.success(request, msg)
    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def recipe_paste(request, uid):
    "Paste recipe stored in the clipboard into project"

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
    attrs.update(stream=recipe.image, name=f"Copy of {recipe.name}", security=recipe.security)
    new_recipe = auth.create_analysis(**attrs)
    new_recipe.save()

    msg = f"Pasted recipe <b>{recipe.name}</b> to project <b>{project.name}</b>."
    msg = mark_safe(msg)
    messages.success(request, msg)
    request.session["recipe_clipboard"] = None

    return redirect(url)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
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

            # Admin users will automatically get authorized.
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
        form = EditCode(user=user, project=project, initial=initial)

    # Bind the JSON to the form.
    recipe = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=dict(name=name))

    # This generates a "fake" unsaved job.
    job = auth.create_job(analysis=analysis, json_data=analysis.json_data, save=False)

    # Create the script for the "fake" job.
    data, script = auth.generate_script(job)

    # Populate the context.
    context = dict(project=project, analysis=analysis, form=form, script=script, recipe=recipe)
    return render(request, 'recipe_code.html', context)


@object_access(type=Project, access=Access.EDIT_ACCESS, url='recipe_list')
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

            recipe = auth.create_analysis(project=project, json_text="{}", template="",
                                          user=request.user, summary=summary, name=name, text=text,
                                          security=security, stream=stream, sticky=sticky)
            recipe.save()
            messages.success(request, "Recipe created")

            return redirect(reverse('recipe_list', kwargs=dict(uid=project.uid)))

    # The url to submit to.
    action_url = reverse('recipe_create', kwargs=dict(uid=project.uid))
    context = dict(project=project, form=form, action_url=action_url, name="New Recipe")

    return render(request, 'recipe_edit.html', context)


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='recipe_view', owner_only=True)
def recipe_edit(request, uid):
    "Edit recipe Info"
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


@object_access(type=Job, access=Access.EDIT_ACCESS, url="job_view", owner_only=True)
def job_edit(request, uid):
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


@object_access(type=Job, access=Access.EDIT_ACCESS, owner_only=True)
def job_delete(request, uid):
    "Change job state to Job.DELETED."

    job = Job.objects.filter(uid=uid).first()
    project = job.project
    job.state = Job.DELETED
    url = reverse('recycle_bin')
    job.save()

    messages.success(request, mark_safe(f"Moved <b>{job.name}</b> to <a href={url}>Recycle Bin</a>."))
    return redirect(reverse("job_list", kwargs=dict(uid=project.uid)))


@object_access(type=Job, access=Access.EDIT_ACCESS, owner_only=True)
def job_restore(request, uid):
    "Change job state from Deleted to Queued."

    job = Job.objects.filter(uid=uid).first()
    job.state = Job.QUEUED
    url = reverse('recycle_bin')
    job.save()

    messages.success(request, mark_safe(f"Moved <b>{job.name}</b> out of <a href={url}>Recycle Bin</a>."))
    return redirect(reverse("job_list", kwargs=dict(uid=job.project.uid)))


@object_access(type=Job, access=Access.READ_ACCESS)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.filter(uid=uid).first()
    project = job.project

    context = dict(job=job, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def block_media_url(request, **kwargs):
    "Block users from accessing media directory using urls"

    messages.error(request, f"Not allowed")
    return redirect(reverse("project_list"))


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_files_list(request, uid, path=''):
    """
    Returns the directory view of the job.
    """
    job = Job.objects.filter(uid=uid).first()
    project = job.project

    form = FileCopyForm(job=job, request=request)
    if request.method == "POST":
        form = FileCopyForm(data=request.POST, job=job, request=request)
        if form.is_valid():
            # Copies files to clipboard
            nfiles = form.save()
            msg = mark_safe(f"Copied <b>{nfiles}</b> file(s) to Clipboard from <b>{job.name}</b>")
            messages.success(request, msg)
            return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))

    context = dict(activate='selection', job=job, project=project, form=form)
    counts = get_counts(project)
    context.update(counts)

    return files_list(request=request, instance=job, path=path,
                      template_name="job_files_list.html", extra_context=context)
