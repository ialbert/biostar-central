import glob
import logging

import mistune
from django.db.models import Q
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse
# from django.template.loader import get_template
from django.utils.safestring import mark_safe
from django.core.paginator import Paginator
from django.views.decorators import csrf, cache
from .const import *
from .decorators import object_access
from .forms import *
from .models import (Project, Data, Analysis, Job, User, Access)

import difflib


def join(*args):
    return os.path.abspath(os.path.join(*args))

# Objects per page when looking at lists
OBJ_PER_PAGE = 10

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


logger = logging.getLogger('engine')


def make_html(text):
    return mistune.markdown(text)

def pages(request, instance):

    paginator = Paginator(instance, OBJ_PER_PAGE)
    page = request.GET.get('page')
    if not page:
        page = 1
    return paginator.page(page)


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
            step = (project.url(), PROJECT_ICON, "Project View", is_active)
        elif icon == DATA_LIST_ICON:
            step = (reverse("data_list", kwargs={'uid': project.uid}), DATA_LIST_ICON, "Data Files", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"File View", is_active)
        elif icon == DATA_UPLOAD:
            step = (reverse("data_view", kwargs={'id': project.id}), DATA_UPLOAD, f"File Upload", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("recipe_list", kwargs={'uid': project.uid}), ANALYSIS_LIST_ICON, "Recipe List", is_active)
        elif icon == ANALYSIS_VIEW_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_VIEW_ICON, "Recipe View", is_active)
        elif icon == ANALYSIS_RUN_ICON:
            step = (reverse("analysis_run", kwargs={'id': analysis.id}), ANALYSIS_RUN_ICON, "Analysis Run", is_active)
        elif icon == ANALYSIS_RECIPE_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_RECIPE_ICON, "Recipe Code", is_active)
        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'uid': project.uid, }), RESULT_LIST_ICON, "Result List", is_active)
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
            step = (reverse("project_view", kwargs={'uid': project.uid}), ADD_USER, "Manage Access", is_active)
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
@csrf.csrf_protect
@cache.never_cache
def project_users(request, uid):
    """
    Manage project users
    """

    project = Project.objects.filter(uid=uid).first()

    # Search query
    q = request.GET.get("q")

    # Users already with access to current project
    users = [access.user for access in project.access_set.all() if access.access > Access.NO_ACCESS]

    # Users that have been searched for.
    targets = []

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ADD_USER],
                               project=project)

    if request.method == "POST":
        form = ChangeUserAccess(data=request.POST)

        if form.is_valid():
            form.change_access()
            messages.success(request, "Changed access to this project")
        else:
            messages.error(request, mark_safe(form.non_field_errors()))

        return redirect(reverse("project_users", kwargs=dict(uid=project.uid)))

    if q:
        targets = User.objects.filter(Q(email__contains=q) | Q(first_name__contains=q))

    current = access_forms(users=users, project=project)
    results = access_forms(users=targets, project=project)

    context = dict(steps=steps, current=current, project=project, results=results)

    return render(request, "project_users.html", context=context)


def project_list(request):

    projects = auth.get_project_list(user=request.user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")
    projects = pages(request, instance=projects)

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])
    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


@object_access(type=Project, access=Access.READ_ACCESS)
@csrf.csrf_protect
@cache.never_cache
def project_view(request, uid):
    user = request.user

    project = Project.objects.filter(uid=uid).first()

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
@csrf.csrf_protect
@cache.never_cache
def project_edit(request, uid):

    project = auth.get_project_list(user=request.user).filter(uid=uid).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON], project=project)

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES, instance=project)
        if form.is_valid():
            form.save()
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

        messages.error(request, mark_safe(form.errors))
        return redirect(reverse("project_edit", kwargs=dict(uid=uid)))

    form = ProjectForm(instance=project)
    context = dict(project=project, steps=steps, form=form)
    return render(request, 'project_edit.html',
                  context)

@csrf.csrf_protect
@cache.never_cache
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
            uid = form.cleaned_data["uid"]
            owner = request.user
            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream, sticky=sticky, privacy=privacy,
                                          uid=uid)
            project.save()
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

        messages.error(request, mark_safe(form.errors))
        return redirect(reverse("project_create"))

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)
    context = dict(steps=steps, form=form)
    return render(request, 'project_create.html', context)


@object_access(type=Project, access=Access.READ_ACCESS)
def data_list(request, uid):
    project = Project.objects.filter(uid=uid).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON],
                               project=project)
    if not project:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    query = Data.objects.filter(project=project).order_by("sticky", "-date")

    data_list = query.all()
    data_count = query.count()
    data_list = pages(request, instance=data_list)

    context = dict(project=project, steps=steps, data_list=data_list, data_count=data_count)
    return render(request, "data_list.html", context)


@object_access(type=Data, access=Access.READ_ACCESS)
@csrf.csrf_protect
@cache.never_cache
def data_view(request, id):

    data = Data.objects.filter(id=id).first()

    if not data:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=data.project, data=data)

    projects = auth.get_project_list(user=request.user)
    projects = projects.exclude(pk=data.project.id).exclude(privacy=Project.PUBLIC)

    # Filter projects by admin access
    cond = Q(access__access__gt=Access.EDIT_ACCESS)
    if request.user.is_authenticated:
        cond = Q(access__user=request.user, access__access__gt=Access.EDIT_ACCESS)
    projects = projects.filter(cond)

    if request.method == "POST":
        form = DataCopyForm(data=request.POST, current=data, request=request)
        name = data.name
        if form.is_valid():
            data = form.save()
            messages.success(request, f"Copied {name} in to {data.project.name}")
        return redirect(reverse("data_view", kwargs=dict(id=data.id)))

    form = DataCopyForm(current=data, request=request)
    context = dict(data=data, steps=steps, projects=projects, form=form)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.EDIT_ACCESS, url='data_view')
@csrf.csrf_protect
@cache.never_cache
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

    form = DataEditForm(instance=data)
    context = dict(data=data, steps=steps, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.UPLOAD_ACCESS, url='data_list')
@csrf.csrf_protect
@cache.never_cache
def data_upload(request, uid):
    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_UPLOAD],
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
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))

        messages.error(request, "Invalid form processing.")
        return redirect(reverse("data_upload", kwargs={'uid': project.uid}))

    form = DataUploadForm()
    context = dict(project=project, steps=steps, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Project, access=Access.READ_ACCESS)
def recipe_list(request, uid):
    """
    Returns the list of analyses for a project id.
    """

    project = Project.objects.filter(uid=uid).first()
    analysis = Analysis.objects.filter(project=project).order_by("-sticky", "-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)
    analysis =  pages(request, instance=analysis)

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "recipe_list.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
@csrf.csrf_protect
@cache.never_cache
def recipe_view(request, id):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON,
                                ANALYSIS_VIEW_ICON],project=analysis.project, analysis=analysis)

    projects = auth.get_project_list(user=request.user)
    projects = projects.exclude(pk=analysis.project.id).exclude(privacy=Project.PUBLIC)

    # Filter projects by admin access
    cond = Q(access__access__gt=Access.EDIT_ACCESS)
    if request.user.is_authenticated:
        cond = Q(access__user=request.user, access__access__gt=Access.EDIT_ACCESS)
    projects = projects.filter(cond)

    if request.method == "POST":
        form = RecipeCopyForm(data=request.POST, analysis=analysis, request=request)
        name = analysis.name
        if form.is_valid():
            analysis = form.save()
            messages.success(request, f"Copied {name} in to {analysis.project.name}")

        return redirect(reverse("recipe_view", kwargs=dict(id=analysis.id)))

    form = RecipeCopyForm(analysis=analysis, request=request)
    context = dict(analysis=analysis, steps=steps, projects=projects, form=form,
                   project=analysis.project)

    return render(request, "recipe_view.html", context)


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

            return redirect(reverse("job_list", kwargs=dict(uid=project.uid)))
        return  redirect(reverse("recipe_run", kwargs=dict(uid=analysis.id)))

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
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

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
            if auth.template_changed(analysis=analysis, template=template):

                # Switch on the untrusted flag when the template changes.
                analysis.security = Analysis.UNDER_REVIEW

                # Set the new template.
                analysis.template = template

            # The SAVE action commits the changes on the analysis.
            if action == 'SAVE':
                analysis.save()
                messages.info(request, "The recipe code has been updated.")
                return redirect(reverse("recipe_view", kwargs=dict(id=analysis.id)))

        return redirect(reverse("recipe_code", kwargs=dict(id=analysis.id)))

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


@object_access(type=Project, access=Access.EDIT_ACCESS, url='recipe_list')
def recipe_create(request, uid):
    """
    Here the id is of the project!
    """
    project = Project.objects.filter(uid=uid).first()

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON], project=project)
    action_url = reverse('recipe_create', kwargs=dict(uid=project.uid))
    back_url = reverse('recipe_list', kwargs=dict(uid=project.uid))

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
        return redirect(action_url)

    form = RecipeForm()
    context = dict(steps=steps, analysis={"name":"New Analysis"},
                   project=project, form=form, action_url=action_url, back_url=back_url)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='recipe_view')
def recipe_edit(request, id):
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

    action_url = reverse('recipe_edit', kwargs=dict(id=analysis.id))
    back_url = reverse('recipe_view', kwargs=dict(id=analysis.id))

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES, instance=analysis)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", kwargs=dict(id=recipe.id)))
        return redirect(action_url)

    form = RecipeForm(instance=analysis)
    context = dict(steps=steps, analysis=analysis, project=project, form=form, action_url=action_url, back_url=back_url)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Project, access=Access.READ_ACCESS, url="project_view")
def job_list(request, uid):
    """
    Returns the list of jobs for a project id.
    """
    project = Project.objects.filter(uid=uid).first()

    if not project:
        messages.error(request, "Jobs not found.")
        # logger.error(f"Jobs for project.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON],
                               project=project)
    jobs = Job.objects.filter(project=project).order_by("-date", "-start_date")

    filter = request.GET.get('filter', '')

    if filter:
        filter = Analysis.objects.filter(id=filter).first()
        jobs = jobs.filter(analysis=filter)

    jobs = pages(request, instance=jobs)

    context = dict(jobs=jobs, steps=steps, project=project, filter=filter)

    return render(request, "job_list.html", context)


@object_access(type=Job, access=Access.EDIT_ACCESS, url="job_view")
def job_edit(request, id):
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON,
                                RESULT_VIEW_ICON], job=job, project=project)

    if request.method == "POST":
        form = JobEditForm(data=request.POST, files=request.FILES, instance=job)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", kwargs=dict(id=job.id)))
        return redirect(reverse("job_edit", kwargs=dict(id=job.id)))

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
        url = settings.MEDIA_URL + job.get_url(path=index)
        return redirect(url)

    return redirect(reverse("job_view", kwargs=dict(id=id)))


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
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
        form = FilesCopyForm(data=request.POST, project=project, job=job)
        if form.is_valid():
            count = form.save()
            messages.success(request, f"Copied {len(count)} file to {project.name}.")
        else:
            messages.warning(request, "Unable to copy files")
        return redirect(reverse("job_view", kwargs=dict(id=job.id)))

    form = FilesCopyForm(project=project)
    context = dict(file_list=file_list, job=job, form=form, steps=steps, project=project, path=path)
    return render(request, "job_files_list.html", context)
