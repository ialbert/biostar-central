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
                     Analysis, Job, make_job)

from engine.const import *
import mistune

def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')


def make_html(text):
    return mistune.markdown(text)


def index(request):

    steps = breadcrumb_builder([HOME_ICON])

    context = dict(steps=steps)

    return render(request, 'index.html', context)


def breadcrumb_builder(icons=[], project=None, analysis=None, data=None, job=None):
    # Works off of icon names

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


@login_required(login_url=LOGIN_URL)
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


@login_required(login_url=LOGIN_URL)
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


@login_required
def data_edit(request, id):

    data = Data.objects.filter(id=id).first()
    project = data.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)

    if request.method == "POST":
        1/0
        form = DataForm(request.POST, instance=data)
        if form.is_valid():
            form.save()

    else:
        1/0
        form = DataForm(instance=data)

    context = dict(data=data, steps=steps, form=form)

    return render(request, 'data_edit.html', context)

@login_required
def data_create(request, id):

    project = Project.objects.filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    if request.method == "POST":

        form = DataForm(request.POST, request.FILES)

        if form.is_valid():

            title = form.cleaned_data["file"]
            text = form.cleaned_data["text"]
            owner = project.owner
            type = form.cleaned_data["type"]

            new_data = Data.objects.create(title=title, text=text,
                                           owner=owner, project=project,
                                           file=request.FILES["file"],
                                           type=type,
                                           )
            new_data.size = f"{os.path.getsize(new_data.file.path)}"
            new_data.save()

            return redirect(reverse("data_list", kwargs={'id':project.id}))

        else:
            form.add_error(None, "Invalid form processing.")
    else:

        form = DataForm()
        context = dict(project=project, steps=steps, form=form)
        return render(request, 'data_create.html', context )


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

    filler = dict(display_type="")

    if spec.get("settings", filler).get("display_type") == "MODEL":
        title = spec["settings"].get("title", analysis.title)
        text = spec["settings"].get("text", analysis.text)
        html = make_html(text)

        return dict(title=title, html=html)
    else:
        return dict()


@login_required
def analysis_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_ICON],
                               project=project, analysis=analysis)
    context = dict()

    if request.method == "POST":
        # get the request.txt and
        # set the mistune
        if request.POST.get("save_or_preview") == "save_to_file":

            form = EditAnalysis(analysis=analysis, data=request.POST)
            if form.is_valid():
                1/0
                return

        elif request.POST.get("save_or_preview") == "preview":

            form = EditAnalysis(analysis=analysis, data=request.POST)

            if form.is_valid():
                form.preview()
                spec = json.loads(form.cleaned_data["text"])
                context = preview_specs(spec, analysis)


        elif request.POST.get("save_or_preview") == "save":

            form = EditAnalysis(analysis=analysis, data=request.POST)

            if form.is_valid():
                form.save()
                spec = json.loads(form.cleaned_data["text"])
                context = preview_specs(spec, analysis)
    else:

        form = EditAnalysis(analysis=analysis)
        spec = json.loads(analysis.json_data)
        context = preview_specs(spec, analysis)

    context.update(dict(project=project, analysis=analysis, steps=steps, form=form))
    print(context["title"])
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


def job_view(request, id):

    # create a directory when clicked.
    job = Job.objects.filter(id=id).first()
    path = job.uid
    url = settings.MEDIA_URL + path+"/"

    return redirect(url)






