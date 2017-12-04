import os, logging, glob
import mistune
from django.conf import settings
# from django.template.loader import get_template
from django.utils.safestring import mark_safe
from django.db.models import Q
from django.contrib import messages
from django.views.decorators import cache
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse
from django import forms

from .forms import *
from .const import *
from .decorators import object_access
from .models import (Project, Data, Analysis, Job, User, Access)

def join(*args):
    return os.path.abspath(os.path.join(*args))

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path =  os.path.abspath(path)
    return path.startswith(__DOCS_DIR)

logger = logging.getLogger('engine')

def make_html(text):
    return mistune.markdown(text)

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
        content = make_html(content)

    title = name.replace("-", " ").replace("_", " ").title()
    context = dict(content=content, title=title, steps=[])
    return render(request, 'info/doc_base.html', context=context)


def index(request):
    steps = breadcrumb_builder([HOME_ICON])
    context = dict(steps=steps)
    return render(request, 'index.html', context)


def breadcrumb_builder(icons=[], project=None, analysis=None, data=None, job=None, user=None):
    """
    This function builds the breadcrumbs on each page.
    """
    if not icons:
        return []

    path = []
    last = icons[-1]
    for icon in icons:
        is_active = icon is last
        if icon == HOME_ICON:
            step = (reverse("index"), HOME_ICON, "Home", is_active)
        elif icon == PROJECT_LIST_ICON:
            step = (reverse("project_list"), PROJECT_LIST_ICON, "Project List", is_active)
        elif icon == PROJECT_ICON:
            step = (reverse("project_view", kwargs={'id': project.id}), PROJECT_ICON, "Project View", is_active)
        elif icon == DATA_LIST_ICON:
            step = (reverse("data_list", kwargs={'id': project.id}), DATA_LIST_ICON, "Data Files", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"File View", is_active)
        elif icon == DATA_UPLOAD:
            step = (reverse("data_view", kwargs={'id': project.id}), DATA_UPLOAD, f"File Upload", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("analysis_list", kwargs={'id': project.id}), ANALYSIS_LIST_ICON, "Recipe List", is_active)
        elif icon == ANALYSIS_VIEW_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_VIEW_ICON, "Recipe View", is_active)
        elif icon == ANALYSIS_RUN_ICON:
            step = (reverse("analysis_run", kwargs={'id': analysis.id}), ANALYSIS_RUN_ICON, "Analysis Run", is_active)
        elif icon == ANALYSIS_RECIPE_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_RECIPE_ICON, "Recipe Code", is_active)
        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'id': project.id, }), RESULT_LIST_ICON, "Result List", is_active)
        elif icon == RESULT_VIEW_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_VIEW_ICON, "Result View", is_active)
        elif icon == USER_ICON:
            step = (reverse("profile"), USER_ICON, f"Profile", is_active)
        elif icon == LOGIN_ICON:
            step = (reverse("login"), LOGIN_ICON, "Login", is_active)
        elif icon == LOGOUT_ICON:
            step = (reverse("login"), LOGOUT_ICON, "Logout", is_active)
        elif icon == INFO_ICON:
            step = (reverse("info"), INFO_ICON, "Information", is_active)
        elif icon == SIGNUP_ICON:
            step = (reverse("signup"), SIGNUP_ICON, "Sign up", is_active)
        elif icon == RESULT_INDEX_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_INDEX_ICON, "Index View", is_active)
        elif icon == ADD_USER:
            step = (reverse("project_view", kwargs={'id': project.id}), ADD_USER, "Manage Access", is_active)
        else:
            continue

        path.append(step)

    return path


@user_passes_test(lambda u: u.is_superuser)
def site_admin(request):
    '''
    Administrative view. Lists the admin project and job.
    '''
    steps = breadcrumb_builder([HOME_ICON])
    projects = Project.objects.all()
    context = dict(steps=steps, projects=projects)
    return render(request, 'admin_index.html', context=context)


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def project_users(request, id):

    project = Project.objects.filter(pk=id).first()
    # Search query, and not_found flag set
    q = not_found = request.GET.get("q")

    # Users already with access to current project
    users = [access.user for access in project.access_set.all() if access.access > Access.NO_ACCESS]
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ADD_USER],
                               project=project)

    if request.method == "POST":
        form = ChangeUserAccess(data=request.POST, project=project, users=users)
        if form.is_valid():
            form.save()
            messages.success(request, "Changed access to this project")
        else:
            # TODO: quick fix for now. not showing up corretly
            messages.error(request, mark_safe(form.non_field_errors()))
        return redirect(reverse("project_users", kwargs=dict(id=id)))

    if q:
        query = User.objects.filter(Q(email__contains=q)|Q(first_name__contains=q))
        # Turn not_found flag off if the query is valid
        not_found = None if query else not_found
        users = query if query else users

    form = ChangeUserAccess(project=project, users=users)
    context = dict(steps=steps, form=form, project=project, not_found=not_found)
    return render(request, "project_users.html", context=context)


def project_list(request):

    projects = auth.get_project_list(user=request.user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])
    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, id):
    user = request.user
    project = Project.objects.filter(id=id).first()

    # Project not found.
    if not project:
        messages.error(request, "Project not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    data_count = Data.objects.filter(project=project).count()
    recipe_count = Analysis.objects.filter(project=project).count()
    result_count = Job.objects.filter(project=project).count()

    if user.is_authenticated():
        access = Access.objects.filter(user=user, project=project).first()
    else:
        access = None

    context = dict(project=project, access=access,
                   data_count=data_count, recipe_count=recipe_count, result_count=result_count,
                   steps=steps)

    return render(request, "project_view.html", context)


@object_access(type=Project, access=Access.EDIT_ACCESS, url='project_view')
def project_edit(request, id):
    project = auth.get_project_list(user=request.user).filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON], project=project)

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES, instance=project)

        if form.is_valid():
            form.save()
        else:
            messages.info(request, f"Invalid form processing")
        return redirect(project.url())

    else:
        form = ProjectForm(instance=project)

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'project_edit.html',
                  context)


def project_create(request):
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    if request.user.is_anonymous:
        messages.warning(request, "You must be logged in to create a project.")
        return redirect(reverse("project_list"))

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

            return redirect(reverse("project_list"))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)
    context = dict(steps=steps, form=form)
    return render(request, 'project_create.html', context)


@object_access(type=Project, access=Access.READ_ACCESS)
def data_list(request, id):

    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON],
                               project=project)
    if not project:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    query = Data.objects.filter(project=project).order_by("sticky", "-date")

    data_list = query.all()
    data_count = query.count()

    context = dict(project=project, steps=steps, data_list=data_list, data_count=data_count)
    return render(request, "data_list.html", context)


# @login_required
@object_access(type=Data, access=Access.READ_ACCESS)
def data_view(request, id):
    data = Data.objects.filter(id=id).first()
    if not data:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=data.project, data=data)
    context = dict(data=data, steps=steps)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.EDIT_ACCESS, url='data_view')
def data_edit(request, id):
    data = Data.objects.filter(id=id).first()
    project = data.project
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)

    if request.method == "POST":
        form = DataEditForm(request.POST, instance=data)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(id=data.id)))
    else:
        form = DataEditForm(instance=data)

    context = dict(data=data, steps=steps, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.UPLOAD_ACCESS, url='data_list')
def data_upload(request, id):
    owner = request.user
    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_UPLOAD ],
                               project=project)

    if request.method == "POST":
        form = DataUploadForm(request.POST, request.FILES)

        if form.is_valid():
            text = form.cleaned_data["text"]
            stream = form.cleaned_data["file"]
            name = stream.name
            auth.create_data(stream=stream, name=name, text=text,
                             user=owner, project=project)
            messages.info(request, "Data upload complete")
            return redirect(reverse("data_list", kwargs={'id': project.id}))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")
        return redirect(reverse("data_upload", kwargs={'id': project.id}))
    else:
        form = DataUploadForm()

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
def analysis_list(request, id):
    """
    Returns the list of analyses for a project id.
    """

    project = Project.objects.filter(id=id).first()
    analyses = Analysis.objects.filter(project=project).order_by("-sticky","-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)
    context = dict(project=project, analyses=analyses, steps=steps)

    return render(request, "analysis_list.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
def recipe_view(request, id):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON],
                               project=project, analysis=analysis)

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "recipe_view.html", context)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view', login_required=True)
def recipe_copy(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    projects = auth.get_project_list(user=request.user)

    # Can't copy into current or public projects
    projects = projects.exclude(pk=analysis.project.id).exclude(privacy=Project.PUBLIC)

    # Filter projects by admin access
    projects = projects.filter(Q(access__user=request.user, access__access__gt=Access.EDIT_ACCESS))
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON],project=analysis.project, analysis=analysis)

    if request.method == "POST":
        form = RecipeCopyForm(data=request.POST, analysis=analysis, user=request.user)
        url = reverse("recipe_copy", kwargs=dict(id=analysis.id))

        if form.is_valid():
            new_analysis = form.save()
            url = reverse("recipe_view", kwargs=dict(id=new_analysis.id))
            messages.success(request, f"Currently in Copy of: {analysis.name}.")

        return redirect(url)

    form = RecipeCopyForm(analysis=analysis, user=request.user)
    context = dict(analysis=analysis, steps=steps, projects=projects, form=form,
                   project=analysis.project, access=Access(access=Access.ADMIN_ACCESS))

    return render(request, "recipe_copy.html", context)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
def recipe_run(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RUN_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, data=request.POST)

        if form.is_valid():

            # The desired name of for the results.
            name = form.cleaned_data.get("name")

            # Generates the JSON data from the bound form field.
            json_data = form.fill_json_data()

            # Create the job from the json.
            job = auth.create_job(analysis=analysis, user=request.user, json_data=json_data, name=name)

            # Spool the job right if UWSGI exists.
            if tasks.HAS_UWSGI:
                jobid = (job.id).to_bytes(5, byteorder='big')
                tasks.execute_job.spool(job_id=jobid)

            return redirect(reverse("job_list", kwargs=dict(id=project.id)))
    else:
        initial = dict(name=analysis.name)
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, steps=steps, form=form)

    return render(request, 'recipe_run.html', context)



@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
def recipe_code(request, id):
    """
    Displays and allows edit on a recipe code.

    Because we allow a preview even for unauthenicated users the view
    is a lot more complicated than a typical DJANO form handler.
    """
    user = request.user

    # There has to be a recipe to work with.
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    name = analysis.name

    # This is the navbat.
    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON],project=project, analysis=analysis)

    if request.method == "POST":
        form = EditCode(user=user, project=project, data=request.POST)

        if form.is_valid():

            # Preview action will let the form cascade through.
            action = form.cleaned_data['action']

            # The changes will commited on SAVE only.
            analysis.json_text = form.cleaned_data['json']

            # We have to check if the template changed.
            template = form.cleaned_data['template']

            # Changes to template will require a review.
            if template != analysis.template:
                # Switch on the untrusted flag when the template changes.
                analysis.security = Analysis.UNDER_REVIEW

                # Set the new template.
                analysis.template = template

            # The SAVE action commits the changes on the analysis.
            if action == 'SAVE':
                analysis.save()
                messages.info(request, "The recipe code has been updated.")
                return redirect(reverse("recipe_view", kwargs=dict(id=analysis.id)))

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
    context = dict(project=project, analysis=analysis, steps=steps, form=form, script=script, recipe=recipe)
    return render(request, 'recipe_code.html', context)


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='analysis_list')
def recipe_create(request, id):
    """
    Here the id is of the project!
    """
    project = Project.objects.filter(id=id).first()

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON], project=project)

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES)
        if form.is_valid():
            recipe = form.save(commit=False)
            # Empty templates may be authorized.
            recipe.security = Analysis.UNDER_REVIEW if recipe.template else Analysis.AUTHORIZED
            recipe.owner = request.user
            recipe.project = project
            recipe.save()
            return redirect(reverse("recipe_view", kwargs=dict(id=recipe.id)))

    form = RecipeForm()
    action_url =reverse('recipe_create', kwargs=dict(id=project.id))
    back_url = reverse('analysis_list', kwargs=dict(id=project.id))
    context = dict(steps=steps, project=project, form=form, action_url=action_url, back_url=back_url)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='recipe_view')
def recipe_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES, instance=analysis, )
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", kwargs=dict(id=analysis.id)))

    form = RecipeForm(instance=analysis)

    action_url = reverse('recipe_edit', kwargs=dict(id=analysis.id))
    back_url = reverse('recipe_view', kwargs=dict(id=analysis.id))

    context = dict(steps=steps, analysis=analysis, project=project, form=form, action_url=action_url, back_url=back_url)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Project, access=Access.READ_ACCESS, url="project_view")
def job_list(request, id):
    """
    Returns the list of jobs for a project id.
    """
    project = Project.objects.filter(id=id).first()

    if not project:
        messages.error(request, "Jobs not found.")
        # logger.error(f"Jobs for project.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON],
                               project=project)
    jobs =Job.objects.filter(project=project).order_by("-date", "-start_date")

    filter = request.GET.get('filter', '')

    if filter:
        filter = Analysis.objects.filter(id=filter).first()
        jobs = jobs.filter(analysis=filter)

    context = dict(jobs=jobs, steps=steps, project=project, filter=filter)

    return render(request, "job_list.html", context)


@object_access(type=Job, access=Access.EDIT_ACCESS, url="job_view")
def job_edit(request, id):
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON,
                                RESULT_VIEW_ICON],job=job, project=project)

    if request.method == "POST":
        form =JobEditForm(data=request.POST, files=request.FILES, instance=job)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", kwargs=dict(id=job.id)))

    form = JobEditForm(instance=job)
    context = dict(steps=steps, job=job, project=project, form=form)

    return render(request, 'job_edit.html', context)


@object_access(type=Job, access=Access.READ_ACCESS)
def job_view(request, id):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON,
                                RESULT_VIEW_ICON], job=job, project=project)

    context = dict(job=job, steps=steps, project=project)
    return render(request, "job_view.html", context=context)


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_result_view(request, id):
    """
    Returns the primary result of a job.
    """
    job = Job.objects.filter(id=id).first()
    index = job.json_data.get("settings", {}).get("index", "")

    if job.state == Job.COMPLETED:

        #TODO:This part is still exposed.
        url = settings.MEDIA_URL + job.get_url(path=index)
        return redirect(url)
    else:
        return redirect(reverse("job_view", kwargs=dict(id=id)))


@object_access(type=Job, access=Access.READ_ACCESS,  url="job_view")
def job_file_view(request, id):
    """
    Returns the directory view of the job.
    """
    job = Job.objects.filter(id=id).first()
    url = settings.MEDIA_URL + job.get_url()

    return redirect(url)


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_files_list(request, id, path=''):
    job = Job.objects.filter(id=id).first()
    project = job.project

    # This is the root of where we can navigate in
    target_path = join(job.path, path)

    if not target_path.startswith(job.path):
        # Attempting to access a file outside of the job directory
        raise Exception(f"target_path {target_path} not in job directory")

    # These are pathlike objects with attributes such as name, is_file
    file_list = list(os.scandir(target_path))

    # Sort by properties
    file_list = sorted(file_list, key=lambda p: (p.is_file(), p.name))

    steps = breadcrumb_builder(
        [PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON, RESULT_VIEW_ICON, RESULT_INDEX_ICON],
        job=job, project=project)

    if request.method == "POST":
        form = DataCopyForm(data=request.POST, project=project, job=job)
        if form.is_valid():
            count = form.save()
            messages.success(request, f"Copied {len(count)} file to {project.name}.")
        else:
            messages.warning(request, "Unable to copy files")
        #TODO: redirection does not make sense really
        return redirect(reverse("job_result_view", kwargs=dict(id=job.id)))

    form = DataCopyForm(project=project)
    context = dict(file_list=file_list, job=job, form=form, steps=steps, project=project, path=path)
    return render(request, "job_files_list.html", context)
