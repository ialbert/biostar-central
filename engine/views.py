import logging
#from django.contrib.auth.decorators import login_required

from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from .settings import BASE_DIR
import os
from.forms import *
from .models import (User, Project, Data,
                     Analysis, Job)

from .util import safe_load


def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')

JSON_SPECFILE =join(BASE_DIR, '..', 'pipeline',
                'templates','metabarcode_qc', 'metabarcode_spec.org.json' )

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
def analysis_list(request):

    analysis = Analysis.objects.order_by("-id")
    if not analysis:
        messages.error(request, "No project found.")
        return redirect("/")

    # True for the active section at the moment
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("analysis_list"), "Analysis List", active)
    ]

    context = dict(analysis=analysis, steps=steps)

    return render(request, "analysis_list.html", context)


#@login_required
def analysis_view(request, id):

    analysis = Analysis.objects.filter(id=id).first()
   # project = analysis.project

    if not analysis:
        messages.error(request, f"Data{id} not found.")

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("analysis_list"),"Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': analysis.id}), f"{analysis.title}", active)
    ]

    context = dict(analysis=analysis, steps=steps)

    return render(request, "analysis_view.html", context)


# also can be seen as results_create
def analysis_run(request, id):

    analysis = Analysis.objects.filter(id=id).first()

    active  = True
    owner = User.objects.all().first()

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("analysis_list"), "Analysis List", not active),
        (reverse("analysis_run", kwargs={'id': analysis.id}), f"{analysis.title}", active),
    ]

    if request.method == "POST":
        # loads the json file into the analysis
        #analysis.load(request.POST.json_file)
        form = RunForm(json_spec=request.POST)
        if form.is_valid:
            form.save()
            filled_json = form.json_spec
            # get the analysis_spec there
            #filled_makefile = form.makefile

            # job makes makefile
            job = Job.objects.get_or_create(json_data= filled_json,
                                            owner=owner,
                                            analysis=analysis)
                                           # makefile=filled_makefile)

            context = dict(job=job, analysis=analysis, steps=steps)
            1/0
            return render(request, "results_list.html", context)

    else:
        #specsfile = open(JSON_SPECFILE)
        # Set json_file in setting.py to avoid loading a file every time a view is served
        analysis.json_spec = JSON_SPECFILE
        analysis.save()
        #specsfile.close()

        form = RunForm(json_spec=json.dumps(safe_load(analysis.json_spec)))
        context = dict(analysis=analysis, steps=steps, form=form)
        return render(request, 'analysis_run.html', context)


def analysis_create(request):

    #analysis = Analysis.objects.filter(id=id).first()

    active  = True
    owner = User.objects.all().first()
    steps = [
        (reverse("index"), "Home", not active),
        (reverse("analysis_list"), "Project List", active),
    ]
    analysis = Analysis.objects.create(owner=owner)
    #analysis.text = analysis.json_spec
    analysis.save()
    if request.method == "POST":

        #json_spec =
        form = EditForm(data=request.POST)
        if form.is_valid():
            # Text field is the json.spec ?
            1/0
            form.save()
    else:

        form = EditForm(json_spec=analysis.json_spec)

    context ={}
    return render(request, 'analysis_create.html', context)


def analysis_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    #analysis.text = analysis.json_spec
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("analysis_list"), "Project List", not active),
        (reverse("analysis_view", kwargs={'id': analysis.id}), f"{analysis.title}", active)
    ]

    if request.method == "POST":

        #form = EditForm(data=request.POST)
        print( request.POST)
        1/0
        if form.is_valid():
            # rewrtite specs file here and save model again.
            1/0
            form.save()

    else:
        #specsfile =
        specs = json.dumps(safe_load(analysis.json_spec), indent=4)
        form = EditForm(json_spec=json.dumps(safe_load(analysis.json_spec)))
        context = dict(analysis=analysis, steps=steps, specs=specs, form=form)
        return render(request, 'analysis_edit.html', context)














