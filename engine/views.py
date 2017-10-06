import logging
#from django.contrib.auth.decorators import login_required
from django.template.loader import get_template
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from .settings import BASE_DIR
import os
from .forms import *
from .models import (User, Project, Data,
                     Analysis, Job)

from .util import make_tmp_jsonfile, rewrite_jsonspecs


def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')

JSON_SPECFILE =join(BASE_DIR, '..', 'pipeline',
                'templates','qc', 'qc_spec.hjson' )

def index(request):

    active = True

    steps = [
        (reverse("index"), "Home", active )
    ]
    context = dict(steps=steps)

    return render(request,'index.html', context)


#@login_required
def project_list(request):

    projects = Project.objects.order_by("-id")

    if not projects:
        messages.error(request, "No project found.")
        return redirect("/")

    # True for the active section at the moment
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", active)
    ]

    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


#@login_required
def project_view(request, id):

    project = Project.objects.filter(id=id).first()

    if not project:
        messages.error(request, f"Project{id} not found.")

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id':project.id}), f"{project.title}", active)
    ]

    context = dict(projects=project, steps=steps)

    return render(request, "project_view.html", context)


def project_edit(request, id):

    project = Project.objects.filter(id=id).first()
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id':project.id}), f"{project.title}", active)
    ]

    if request.method == "POST":

        form = ProjectForm(request.POST, instance=project)
        if form.is_valid():
            form.save()

    else:
        form = ProjectForm(instance=project)

    context = dict(projects=project, steps=steps, form=form)
    return render(request, 'project_edit.html',
                      context)


def project_create(request):

    active = True
    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_create"), "Create Project", active)
    ]
    if request.method == "POST":
        # create new projects here ( just populate metadata ).

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
    if not project:
        messages.error(request, "No data found for this project.")

    active = True
    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", not active),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", active),
    ]

    context = dict(project=project, steps=steps)

    return render(request, "data_list.html", context)


#@login_required
def data_view(request, id):

    data = Data.objects.filter(id=id).first()
    project = data.project

    if not data:
        messages.error(request, f"Data{id} not found.")

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", not active),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", not active),
        (reverse("data_view", kwargs={'id': data.id}), f"{data.title}", active)
    ]

    context = dict(data=data, steps=steps)

    return render(request, "data_view.html", context)


def data_edit(request, id):

    data = Data.objects.filter(id=id).first()
    project = data.project

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", not active),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", not active),
        (reverse("data_view", kwargs={'id': data.id}), f"{data.title}", active)
    ]

    if request.method == "POST":

        form = DataForm(request.POST, instance=data)
        if form.is_valid():
            form.save()

    else:
        form = DataForm(instance=data)

    context = dict(data=data, steps=steps, form=form)

    return render(request, 'data_edit.html', context)


def data_create(request, id):

    project = Project.objects.filter(id=id).first()
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", not active),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", not active),
        (reverse("data_create", kwargs={'id': project.id}), "Create Data", active)
    ]

    if request.method == "POST":

        form = DataForm(data=request.POST)

        if form.is_valid():

            title = form.cleaned_data["title"]
            text = form.cleaned_data["text"]
            owner = User.objects.all().first()
            new_data = Data.objects.create(title=title, text=text,
                                           owner=owner, project=project)
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

    analysis = Analysis.objects.order_by("-id")
    project = Project.objects.filter(id=id).first()

    if not analysis:
        messages.error(request, "No Analysis found.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    # True for the active section at the moment
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id':project.id}), "Analysis List", active)
    ]

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_list.html", context)


#@login_required
def analysis_view(request, id, id2):

    analysis = Analysis.objects.filter(id=id2).first()
    project = Project.objects.filter(id=id).first()
   # project = analysis.project

    if not analysis:
        messages.error(request, "Analysis not found.")
        return redirect(reverse("analysis_list", kwargs={'id':project.id}))

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id':project.id}), "Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': project.id, 'id2':analysis.id}),
         f"{analysis.title}", active)
    ]

    context = dict(project=project, analysis=analysis, steps=steps)

    return render(request, "analysis_view.html", context)


def analysis_run(request, id, id2):

    project = Project.objects.filter(id=id).first()
    analysis = Analysis.objects.filter(id=id2).first()
    owner = User.objects.all().first()

    active  = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id': project.id}), "Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': project.id, 'id2': analysis.id}),
         f"{analysis.title}", active)
    ]

    if request.method == "POST":

        steps = [
            (reverse("index"), "Home", not active),
            (reverse("project_list"), "Project List", not active),
            (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
            (reverse("jobs_list", kwargs={'id': project.id}), "Results List", active),
        ]

        form = RunForm(data=request.POST, json_spec=analysis.json_spec)

        if form.is_valid():

            # Fill "value" with what the user picked.
            filled_json = safe_load(analysis.json_spec)
            # add template to spec
            for field in filled_json:
                if field in form.cleaned_data:
                    # Mutates field value in spec
                    data = filled_json[field]
                    data["value"] = form.cleaned_data[field]

            # path comes from spec
            makefile_template = get_template("qc/qc_makefile.html").template.source

            job = Job.objects.get_or_create(json_data=filled_json,
                                            owner=owner,
                                            analysis=analysis,
                                            project=project,
                                            makefile_template=makefile_template)

            context = dict(jobs=project.job_set, job=job, steps=steps)
            # return redirect(reverse(jobs_list))
            return render(request, "jobs_list.html", context)

    else:

        analysis.json_spec = JSON_SPECFILE
        analysis.save()
        # RunAnalysis(analysis=analysis)
        form = RunForm(json_spec=analysis.json_spec)
        context = dict(project=project, analysis=analysis, steps=steps, form=form)
        return render(request, 'analysis_run.html', context)


def analysis_edit(request, id, id2):

    analysis = Analysis.objects.filter(id=id2).first()
    project = Project.objects.filter(id=id).first()

    analysis.json_spec = JSON_SPECFILE
    # should only save in analysis_create!
    analysis.save()
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id': project.id}), "Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': project.id, 'id2': analysis.id}),
         f"{analysis.title}", active)
    ]

    if request.method == "POST":

        if request.POST.get("save_or_preview") == "save":
            specs = analysis.json_spec

            rewrite_jsonspecs(request.POST.get("text"), specs)
            form = EditForm(json_spec=specs)

        elif request.POST.get("save_or_preview") == "preview":

            tmp_specs = make_tmp_jsonfile(request.POST.get("text"), analysis.id)
            form = EditForm(json_spec=tmp_specs)

    else:
        form = EditForm(json_spec=analysis.json_spec)

    context = dict(project=project, analysis=analysis, steps=steps, form=form)
    return render(request, 'analysis_edit.html', context)


def jobs_list(request, id):

    project = Project.objects.filter(id=id).first()

    if not project.job_set:
        messages.error(request, "No Jobs found.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("jobs_list", kwargs={'id': project.id}), "Results List", active),
    ]

    context = dict(jobs=project.job_set, steps=steps)

    return render(request, "jobs_list.html", context)


def job_view(request):
    return






