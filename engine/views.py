import logging
#from django.contrib.auth.decorators import login_required
from django.template.loader import get_template
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
#from django.template import Context, loader
from django.conf import settings
from django.conf.urls.static import serve
import os
from .forms import *
from .models import (User, Project, Data,
                     Analysis, Job)

from . import util


def join(*args):
    return os.path.abspath(os.path.join(*args))


logger = logging.getLogger('engine')

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

    if not projects.all():
        messages.error(request, "No projects found.")
        return redirect(reverse("index"))

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
        messages.error(request, "Project not found.")

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
        messages.error(request, "Data not found.")
        return redirect(reverse("data_view", kwargs={'id': data.id}))

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

        form = DataForm(request.POST, request.FILES)

        if form.is_valid():

            title = form.cleaned_data["title"]
            text = form.cleaned_data["text"]
            owner = User.objects.all().first()
            type = form.cleaned_data["type"]

            # maybe specifiy path here instead of in save method
            new_data = Data.objects.create(title=title, text=text,
                                           owner=owner, project=project,
                                           file=request.FILES["file"],
                                           type=type)
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
    analysis.save()
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

        form = RunAnalysis(data=request.POST, analysis=analysis.spec_source)

        if form.is_valid():
            filled_json = util.safe_loads(analysis.spec_source)

            for field in filled_json:
                data = filled_json[field]
                if field in form.cleaned_data:
                    # Mutates field value in spec_source
                    data["value"] = form.cleaned_data[field]

            template_path = filled_json["template"]["value"]

            title = form.cleaned_data["title"]
            if form.cleaned_data["title"] == "Title":
                title = analysis.title

            makefile_template = get_template(template_path).template.source

            job = Job.objects.create(json_data=filled_json,
                                            owner=owner,
                                            analysis=analysis,
                                            project=project,
                                            makefile_template=makefile_template,
                                            title=title)

            job.save()
            return redirect(reverse("jobs_list", kwargs=dict(id=project.id)))

    else:
        form = RunAnalysis(analysis=analysis.spec_source)
        context = dict(project=project, analysis=analysis, steps=steps, form=form)
        return render(request, 'analysis_run.html', context)


def analysis_edit(request, id, id2):

    analysis = Analysis.objects.filter(id=id2).first()
    project = Project.objects.filter(id=id).first()

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

        if request.POST.get("save_or_preview") == "save_to_file":

            # NEED TO VALIDATE BEFORE OVERRIDING FILE
            spec_file = analysis.spec_origin
            util.rewrite_specs(request.POST.get("text"), spec_file)
            # Save the rewriteen spec_file into spec_origin field
            analysis.save()

            form = EditAnalysis(analysis=analysis.spec_origin)

        elif request.POST.get("save_or_preview") == "preview":

            form = EditAnalysis(analysis=request.POST.get("text"))

        elif request.POST.get("save_or_preview") == "save":

            form = EditAnalysis(analysis=request.POST.get("text"))
            spec = util.safe_loads(analysis.spec_source)
            filler = dict(display_type='')

            if spec.get("analysis_spec", filler)["display_type"] == "MODEL":
                analysis.title = spec["analysis_spec"].get("title", analysis.title)
                analysis.text = spec["analysis_spec"]["value"]

            analysis.save(spec_source=request.POST.get("text"))
    else:

        form = EditAnalysis(analysis=analysis.spec_source)

    context = dict(project=project, analysis=analysis, steps=steps, form=form)
    return render(request, 'analysis_edit.html', context)


def jobs_list(request, id):

    project = Project.objects.filter(id=id).first()

    if not project.job_set.all():
        messages.error(request, "No jobs found for this project.")
        return redirect(reverse("project_view", kwargs={'id': project.id}))

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("jobs_list", kwargs={'id': project.id}), "Results List", active),
    ]

    jobs = project.job_set.order_by("-id")

    context = dict(jobs=jobs, steps=steps)

    return render(request, "jobs_list.html", context)


def media_files(request, path):

    return serve(request, path=path, document_root=settings.MEDIA_DIR, show_indexes=True)
def job_view(request, id):

    # create a directory when clicked.
    job = Job.objects.filter(id=id).first()
    path = job.uid
    url = settings.MEDIA_URL + path+"/"
    dir_root = join(job.path, "..", "..")

    #return serve(request, path=url, document_root=dir_root, show_indexes=True)

    print(url)
    return redirect(url)






