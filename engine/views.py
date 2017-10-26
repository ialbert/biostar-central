# import os

import mistune
from django.conf import settings
# from django.template.loader import get_template
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.urls import reverse

from . import tasks
from .forms import *
from .models import (User, Project, Data,
                     Analysis, Job, get_datatype)


def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')


def make_html(text):
    return mistune.markdown(text)


def info(request):
    steps = breadcrumb_builder([HOME_ICON, INFO_ICON])
    context = dict(steps=steps, info=make_html(INFO))
    return render(request, 'info.html', context=context)


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
            step = (reverse("data_list", kwargs={'id': project.id}), DATA_LIST_ICON, "Data List", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"Data View", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("analysis_list", kwargs={'id': project.id}), ANALYSIS_LIST_ICON, "Analysis List", is_active)
        elif icon == ANALYSIS_ICON:
            step = (reverse("analysis_view", kwargs={'id': analysis.id}), ANALYSIS_ICON, "Analysis View", is_active)
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


# TODO: Replace with the manager thing.
def project_list(request):
    projects = Project.objects.order_by("-id")
    if not request.user.is_superuser:
        projects = projects.filter(type=Project.USER).all()

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
        form = ProjectForm(request.POST, instance=project)
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
        form = ProjectForm(data=request.POST)
        if form.is_valid():
            name = form.cleaned_data["name"]
            text = form.cleaned_data["text"]
            owner = User.objects.all().first()
            project = Project.objects.create(name=name, text=text, owner=owner)
            project.save()
            return redirect(reverse("project_list"))
        else:
            form.add_error(None, "Invalid form processing.")
    else:
        initial = dict(name="Project Name", text="project description", summary="project summary")
        form = ProjectForm(initial=initial)
        context = dict(steps=steps, form=form)
        return render(request, 'project_create.html',
                      context)


# @login_required
def data_list(request, id):
    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON],
                               project=project)
    context = dict(project=project, steps=steps)
    return render(request, "data_list.html", context)


# @login_required
def data_view(request, id):
    data = Data.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=data.project, data=data)
    context = dict(data=data, steps=steps)

    return render(request, "data_view.html", context)


def remove_file(file):
    try:
        os.remove(file.path)
    except FileNotFoundError:
        pass
    return


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
            data_type = get_datatype(name)
            project.create_data(stream=stream, name=name, data_type=data_type, text=text,
                                owner=owner)

            return redirect(reverse("data_list", kwargs={'id': project.id}))

        else:
            form.add_error(None, "Invalid form processing.")
    else:
        form = DataUploadForm()
        context = dict(project=project, steps=steps, form=form)
        return render(request, 'data_upload.html', context)


# @login_required
def analysis_list(request, id):
    """
    Returns the list of analyses for a project id.
    """
    # filter according to user.

    project = Project.objects.filter(id=id).first()
    analysis = Analysis.objects.filter(project=project).order_by("-id")

    if not request.user.is_superuser:
        analysis = analysis.filter(type=Analysis.USER).all()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)
    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_list.html", context)


def analysis_view(request, id):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_ICON],
                               project=project, analysis=analysis)
    if request.method == "POST":
        form = ExportAnalysis(data=request.POST, analysis=analysis)
        if form.is_valid():
            project, analysis = form.export()
            return redirect(reverse("analysis_list", kwargs={"id": project.id}))
    else:
        form = ExportAnalysis(analysis=analysis)
    context = dict(project=project, analysis=analysis, steps=steps,
                   form=form)

    return render(request, "analysis_view.html", context)


@login_required
def analysis_run(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":
        form = RunAnalysis(data=request.POST, analysis=analysis)

        if form.is_valid():
            name = form.cleaned_data.get("name")
            filled_json = form.process()
            json_text = hjson.dumps(filled_json)
            job = analysis.create_job(owner=analysis.owner, json_text=json_text, name=name,
                                      type=analysis.type)
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
    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_ICON],
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
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON],
                               project=project)

    jobs = project.job_set.order_by("-id")
    if not request.user.is_superuser:
        jobs = jobs.filter(type=Job.USER).all()

    context = dict(jobs=jobs, steps=steps, project=project)

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

    if job.state == Job.FINISHED:
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

        form = DataCopyForm(data=request.POST, project=project)
        if form.is_valid():
            count = form.process()
            messages.success(request, f"Copied {count} file to {project.name}.")
        else:
            messages.warning(request, "Unable to copy files")

    else:
        form = DataCopyForm(project=project)

    context = dict(file_list=file_list, job=job, form=form, steps=steps, project=project, path=path)
    return render(request, "job_files_list.html", context)
