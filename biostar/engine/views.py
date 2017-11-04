# import os

import mistune
from django.conf import settings
# from django.template.loader import get_template
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse

from . import auth
from . import tasks
from .const import *
from .forms import *
from .models import (User, Project, Data,
                     Analysis, Job)


def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')


def make_html(text):
    return mistune.markdown(text)


def info(request):
    tmp = "Store and analyze metagenomic data"
    steps = breadcrumb_builder([HOME_ICON, INFO_ICON])
    context = dict(steps=steps, info=make_html(tmp))
    return render(request, 'docs/info.html', context=context)


def index(request):
    steps = breadcrumb_builder([HOME_ICON])
    context = dict(steps=steps)
    return render(request, 'index.html', context)


def breadcrumb_builder(icons=[], project=None, analysis=None, data=None, job=None, user=None):
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
            step = (reverse("data_list", kwargs={'id': project.id}), DATA_LIST_ICON, "File List", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"File View", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("analysis_list", kwargs={'id': project.id}), ANALYSIS_LIST_ICON, "Analysis List", is_active)
        elif icon == ANALYSIS_VIEW_ICON:
            step = (
            reverse("analysis_view", kwargs={'id': analysis.id}), ANALYSIS_VIEW_ICON, "Analysis View", is_active)
        elif icon == ANALYSIS_RUN_ICON:
            step = (reverse("analysis_run", kwargs={'id': analysis.id}), ANALYSIS_RUN_ICON, "Analysis Run", is_active)
        elif icon == ANALYSIS_RECIPE_ICON:
            step = (reverse("analysis_recipe", kwargs={'id': analysis.id}), ANALYSIS_RECIPE_ICON, "Analysis Recipe", is_active)

        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'id': project.id, }), RESULT_LIST_ICON, "Result List", is_active)
        elif icon == RESULT_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_ICON, "Job View", is_active)
        elif icon == RESULT_VIEW_ICON:
            step = (reverse("job_detail_view", kwargs={'id': job.id, }), RESULT_ICON, "Job Status", is_active)
        elif icon == USER_ICON:
            step = (reverse("profile", kwargs={'id': user.id, }), USER_ICON, f"Profile", is_active)
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


def project_list(request):
    projects = Project.objects.order_by("-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


# @login_required
def project_view(request, id):
    project = Project.objects.filter(id=id).first()

    # Project not found.
    if not project:
        messages.error(request, "Project not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    context = dict(project=project, steps=steps)

    return render(request, "project_view.html", context)


@login_required
def project_edit(request, id):
    project = Project.objects.filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON], project=project)

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES, instance=project)

        if form.is_valid():
            form.save()
            return redirect(project.url())

    else:
        form = ProjectForm(instance=project)

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'project_edit.html',
                  context)


@login_required
def project_create(request):
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    if request.method == "POST":
        # create new projects here ( just populates metadata ).
        form = ProjectForm(request.POST, request.FILES)
        if form.is_valid():

            name = form.cleaned_data["name"]
            text = form.cleaned_data["text"]
            summary = form.cleaned_data["summary"]
            stream = form.cleaned_data["image"]
            owner = request.user

            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream)
            project.save()

            return redirect(reverse("project_list"))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)
    context = dict(steps=steps, form=form)
    return render(request, 'project_create.html',
                  context)


# @login_required
def data_list(request, id):
    project = Project.objects.get(id=id)
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON],
                               project=project)
    if not project:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    data_list = Data.objects.filter(project=project).order_by("-date")
    context = dict(project=project, steps=steps, data_list=data_list)
    return render(request, "data_list.html", context)


# @login_required
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


@login_required
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


@login_required
def data_upload(request, id):
    owner = request.user

    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    if request.method == "POST":
        form = DataUploadForm(request.POST, request.FILES)

        if form.is_valid():
            text = form.cleaned_data["text"]
            stream = form.cleaned_data["file"]
            name = stream.name
            data_type = form.cleaned_data["data_type"]

            auth.create_data(stream=stream, name=name, data_type=data_type, text=text,
                             user=owner, project=project)
            messages.info(request, "Data upload complete")
            return redirect(reverse("data_list", kwargs={'id': project.id}))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")

    else:
        form = DataUploadForm()

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'data_upload.html', context)


def analysis_list(request, id):
    """
    Returns the list of analyses for a project id.
    """
    # filter according to user.

    project = Project.objects.filter(id=id).first()
    analyses = Analysis.objects.filter(project=project).order_by("-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)
    context = dict(project=project, analyses=analyses, steps=steps)

    return render(request, "analysis_list.html", context)


def analysis_view(request, id):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON],
                               project=project, analysis=analysis)

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_view.html", context)


def analysis_recipe(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RECIPE_ICON],
                               project=analysis.project, analysis=analysis)

    context=dict(analysis=analysis, steps=steps)
    return render(request, "analysis_recipe.html", context)


def analysis_copy(request, id):

    #TODO: will use a factory.py function for generating projects field when adding new features
    analysis = Analysis.objects.filter(id=id).first()
    projects = Project.objects.all()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RECIPE_ICON],
                               project=analysis.project, analysis=analysis)

    if request.method == "POST":

        form = AnalysisCopyForm(data=request.POST, analysis=analysis)
        if form.is_valid():
            count = form.process()
            messages.success(request, f"Copied current analysis to {count} project(s).")
    else:
        form = AnalysisCopyForm(analysis=analysis)

    context=dict(analysis=analysis, steps=steps, projects=projects, form=form)
    return render(request, "analysis_copy.html", context)


def analysis_run(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    project = analysis.project

    steps = breadcrumb_builder( [HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
         ANALYSIS_VIEW_ICON, ANALYSIS_RUN_ICON],
        project=project, analysis=analysis)

    if request.method == "POST":
        form = RunAnalysis(data=request.POST, analysis=analysis)

        if request.user.is_anonymous():
            form.add_error(None, "You must be logged in to run an analysis")

        if form.is_valid():
            name = form.cleaned_data.get("name")
            filled_json = form.process()
            json_text = hjson.dumps(filled_json)
            job = auth.create_job(analysis=analysis, user=analysis.owner, json_text=json_text, name=name,
                                  )
            logger.info(tasks.HAS_UWSGI)
            if tasks.HAS_UWSGI:
                jobid = (job.id).to_bytes(5, byteorder='big')
                tasks.execute_job.spool(job_id=jobid)

            return redirect(reverse("job_list", kwargs=dict(id=project.id)))
    else:
        initial = dict(name=analysis.name)
        form = RunAnalysis(analysis=analysis, initial=initial)

    context = dict(project=project, analysis=analysis, steps=steps, form=form)
    return render(request, 'analysis_run.html', context)


def preview_specs(spec, analysis):
    if spec.get("settings"):
        name = spec["settings"].get("name", analysis.name)
        help = spec["settings"].get("help", analysis.text)
        # summary = spec["settings"].get("summary", analysis.text)
        html = make_html(help)

        return dict(name=name, html=html)
    else:
        return dict()


def process_analysis_edit(method, analysis, form):
    form_method_map = {'preview': form.preview,
                       'save': form.save}
    spec = dict()
    if form.is_valid():
        form_method_map[method]()
        spec = hjson.loads(form.cleaned_data["text"])

    return preview_specs(spec, analysis)


@login_required
def analysis_edit(request, id):
    analysis = Analysis.objects.filter(id=id).first()
    # filter according to user
    project = analysis.project
    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":
        form = EditAnalysisForm(analysis=analysis, data=request.POST)
        method = request.POST.get("save_or_preview")
        context = process_analysis_edit(method, analysis, form)

    else:
        form = EditAnalysisForm(analysis=analysis)
        spec = hjson.loads(analysis.json_text)
        context = preview_specs(spec, analysis)

    context.update(dict(project=project, analysis=analysis, steps=steps, form=form))
    return render(request, 'analysis_edit.html', context)


def job_list(request, id):
    """
    Returns the list of jobs for a project id.
    """
    # filter according to type
    project = Project.objects.filter(id=id).first()

    if not project:
        messages.error(request, "Jobs not found.")
        # logger.error(f"Jobs for project.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON],
                               project=project)

    jobs = project.job_set.order_by("-id")

    filter = request.GET.get('filter', '')

    if filter:
        filter = Analysis.objects.filter(id=filter).first()
        jobs = jobs.filter(analysis=filter)

    context = dict(jobs=jobs, steps=steps, project=project, filter=filter)

    return render(request, "job_list.html", context)


def job_view(request, id):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON, RESULT_VIEW_ICON],
                               job=job, project=project)

    context = dict(job=job, steps=steps)
    return render(request, "job_view.html", context=context)


def job_result_view(request, id):
    """
    Returns the primary result of a job.
    """
    job = Job.objects.filter(id=id).first()
    index = job.json_data.get("settings", {}).get("index", "")

    if job.state == Job.COMPLETED:
        url = settings.MEDIA_URL + job.get_url(path=index)
        return redirect(url)
    else:
        return job_view(request, job.id)


def job_file_view(request, id):
    """
    Returns the directory view of the job.
    """
    job = Job.objects.filter(id=id).first()
    url = settings.MEDIA_URL + job.get_url()

    return redirect(url)


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

    steps = breadcrumb_builder(
        [PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON, RESULT_VIEW_ICON, RESULT_INDEX_ICON],
        job=job, project=project)

    if request.method == "POST":

        form = DataCopyForm(data=request.POST, project=project, job=job)
        if form.is_valid():
            count = form.process()
            messages.success(request, f"Copied {count} file to {project.name}.")
        else:
            messages.warning(request, "Unable to copy files")

    else:
        form = DataCopyForm(project=project)

    context = dict(file_list=file_list, job=job, form=form, steps=steps, project=project, path=path)
    return render(request, "job_files_list.html", context)
