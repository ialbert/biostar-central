import os
import hjson as json
import logging
from django.contrib.auth.decorators import login_required
from django.template.loader import get_template
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from django.conf import settings
from .forms import *
from .models import (User, Project, Data,
                     Analysis, Job, make_job, get_datatype)
from django.core.files import File

from engine.const import *
import mistune


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


def breadcrumb_builder(icons=[], project=None, analysis=None, data=None, job=None):

    if not icons:
        return []

    path = []
    last = icons[-1]
    for icon in icons:
        is_active = icon is last
        if icon == HOME_ICON:
            step = (reverse("index"), HOME_ICON, "Home", is_active )
        elif icon == PROJECT_LIST_ICON:
            step = (reverse("project_list"), PROJECT_LIST_ICON, "Project List", is_active)
        elif icon == PROJECT_ICON:
            step = (reverse("project_view", kwargs={'id': project.id}), PROJECT_ICON, f"{project.title}", is_active )
        elif icon == DATA_LIST_ICON:
            step = (reverse("data_list", kwargs={'id': project.id}), DATA_LIST_ICON, "Data List", is_active )
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"{data.title}", is_active )
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("analysis_list", kwargs={'id': project.id}), ANALYSIS_LIST_ICON, "Analysis List", is_active )
        elif icon == ANALYSIS_ICON:
            step = (reverse("analysis_view", kwargs={'id': analysis.id}), ANALYSIS_ICON, f"{analysis.title}", is_active )
        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'id': project.id,}), RESULT_LIST_ICON, "Result List",is_active)
        elif icon == RESULT_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_ICON, f"{job.title}", is_active)
        elif icon == LOGIN_ICON:
            step = (reverse("login"), LOGIN_ICON, "Login", is_active)
        elif icon == LOGOUT_ICON:
            step = (reverse("login"), LOGOUT_ICON, "Logout", is_active)
        elif icon == INFO_ICON:
            step = (reverse("info"), INFO_ICON, "Information", is_active)
        elif icon == SIGNUP_ICON:
            step = (reverse("signup"), SIGNUP_ICON, "Sign up", is_active)
        elif icon == RESULT_VIEW_ICON:
            step = (reverse("job_detail_view", kwargs={'id': job.id,}), RESULT_ICON, f"{job.title}", is_active)
        else:
            continue

        path.append(step)

    return path


#@login_required
def project_list(request):

    projects = Project.objects.order_by("-id")

    if not projects.all():
        messages.error(request, "No projects found.")
        return redirect(reverse("index"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


#@login_required
def project_view(request, id):

    project = Project.objects.filter(id=id).first()

    if not project:
        messages.error(request, "Project not found.")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON], project=project)

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

    else:
        form = ProjectForm(instance=project)

    context = dict(projects=project, steps=steps, form=form)
    return render(request, 'project_edit.html',
                      context)


@login_required
def project_create(request):

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    if request.method == "POST":
        # create new projects here ( just populates metadata ).

        form = ProjectForm(data=request.POST)

        if form.is_valid():

            title = form.cleaned_data["title"]
            text = form.cleaned_data["text"]
            owner = User.objects.all().first()
            project = Project.objects.create(title=title, text=text, owner=owner)
            project.save()

            return redirect(reverse("project_list"))
        else:
            form.add_error(None, "Invalid form processing.")
    else:

        form = ProjectForm()
        context = dict(steps=steps, form=form)
        return render(request, 'project_create.html',
                      context)


#@login_required
def data_list(request, id):

    project = Project.objects.filter(id=id).first()

    if not project.data_set.all():
        messages.error(request, "No data found for this project.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    steps = breadcrumb_builder([HOME_ICON,  PROJECT_ICON, DATA_LIST_ICON], project=project)

    context = dict(project=project, steps=steps)

    return render(request, "data_list.html", context)


#@login_required
def data_view(request, id):

    data = Data.objects.filter(id=id).first()
    project = data.project

    if not data:
        messages.error(request, "Data not found.")
        return redirect(reverse("data_view", kwargs={'id': data.id}))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)

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
    initial = dict(text=data.text)

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)

    if request.method == "POST":

        form = DataEditForm(request.POST, initial=initial)

        if form.is_valid():

            data.text = form.cleaned_data["text"]

            data.save()

    else:
        form = DataEditForm(initial=initial)

    context = dict(data=data, steps=steps, form=form)

    return render(request, 'data_edit.html', context)


@login_required
def data_upload(request, id):

    project = Project.objects.filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    if request.method == "POST":

        form = DataUploadForm(request.POST, request.FILES)

        if form.is_valid():

            title = form.cleaned_data["file"]
            data_file = str(form.cleaned_data["file"])
            file = File(form.cleaned_data["file"])
            owner = User.objects.filter(email=request.user).first() or project.owner
            text = form.cleaned_data["text"]
            data_type = get_datatype(file)
            size = f"{file.size}"
            data = Data(title=title, owner=owner, text=text, project=project, data_type=data_type, size=size)
            data.file.save(data_file, file, save=True)
            data.save()

            return redirect(reverse("data_list", kwargs={'id':project.id}))

        else:
            form.add_error(None, "Invalid form processing.")
    else:
        form = DataUploadForm()
        context = dict(project=project, steps=steps, form=form)
        return render(request, 'data_upload.html', context )

#@login_required
def analysis_list(request, id):

    project = Project.objects.filter(id=id).first()
    analysis = Analysis.objects.filter(project=project).order_by("-id")

    if not analysis:
        messages.error(request, "No Analysis found.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    steps = breadcrumb_builder([HOME_ICON,  PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_list.html", context)


def analysis_view(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    if not analysis:
        messages.error(request, "Analysis not found.")
        return redirect(reverse("analysis_list", kwargs={'id':project.id}))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON, ANALYSIS_ICON],
                               project=project, analysis=analysis)

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_view.html", context)


@login_required
def analysis_run(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    owner = analysis.owner
    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON,  ANALYSIS_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":

        form = RunAnalysis(data=request.POST, analysis=analysis)

        if form.is_valid():

            filled_json = form.process()
            json_data = json.dumps(filled_json)
            title = form.cleaned_data.get("title")
            job = make_job(owner=owner, analysis=analysis, project=project,
                           json_data=json_data, title=title)

            job.save()

            return redirect(reverse("job_list", kwargs=dict(id=project.id)))

    else:
        initial = dict(title=analysis.title)
        form = RunAnalysis(analysis=analysis, initial=initial)
        context = dict(project=project, analysis=analysis, steps=steps, form=form)
        return render(request, 'analysis_run.html', context)


def preview_specs(spec, analysis):

    if spec.get("settings"):
        title = spec["settings"].get("title", analysis.title)
        help = spec["settings"].get("help", analysis.text)
        #summary = spec["settings"].get("summary", analysis.text)
        html = make_html(help)

        return dict(title=title, html=html)
    else:
        return dict()


def process_analysis_edit(method, analysis, form):

    form_method_map = {'preview':form.preview,
                       'save':form.save,
                       'save_to_file': form.save_to_file}
    spec = dict()
    if form.is_valid():

        form_method_map[method]()
        spec = json.loads(form.cleaned_data["text"])

    return preview_specs(spec, analysis)


@login_required
def analysis_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":

        form = EditAnalysisForm(analysis=analysis, data=request.POST)
        method = request.POST.get("save_or_preview")
        context = process_analysis_edit(method, analysis, form)

    else:

        form = EditAnalysisForm(analysis=analysis)
        spec = json.loads(analysis.json_data)
        context = preview_specs(spec, analysis)

    context.update(dict(project=project, analysis=analysis, steps=steps, form=form))
    return render(request, 'analysis_edit.html', context)


def jobs_list(request, id):

    project = Project.objects.filter(id=id).first()

    if not project.job_set.all():
        messages.error(request, "No jobs found for this project.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON ],
                               project=project)

    jobs = project.job_set.order_by("-id")

    context = dict(jobs=jobs, steps=steps, project=project)

    return render(request, "jobs_list.html", context)


def media_index(request):

    context = dict()

    return render(request, "media_index.html", context)

@login_required
def job_view(request, id):

    # create a directory when clicked.
    job = Job.objects.filter(id=id).first()
    path = f"JOB{job.uid}"
    url = settings.MEDIA_URL + path+"/"

    return redirect(url)

def job_detail_view(request, id):

    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON, RESULT_LIST_ICON, RESULT_VIEW_ICON ],
                               job=job, project=project)

    context = dict(job=job, steps=steps)
    return render(request, "job_view.html", context=context)




