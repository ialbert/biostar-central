import uuid
import logging

from django.contrib.auth.decorators import login_required
from django.contrib.auth import login, authenticate, logout
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse
from ratelimit.decorators import ratelimit


from.forms import *
from .models import (User, Project, Data,
                     Analysis, Result, Job)


logger = logging.getLogger('engine')


def index(request):

    active = True


    steps = [
        (reverse("index"), "Home", active )
    ]
    can_create_project = True
    context = dict(steps=steps, can_create_project=can_create_project)

    return render(request,'index.html', context)


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_signup(request):

    if request.method == 'POST':
        form = SignUpForm(request.POST)

        if form.is_valid():

            email = form.cleaned_data.get('email')
            password = form.cleaned_data.get('password1')
            name = email.split("@")[0]
            
            user = User.objects.create(username=get_uuid(), email=email, 
                                       first_name=name)
            user.set_password(password)
            user.save()

            login(request, user)
            logger.info(f"Signed up and logged in user.id={user.id}, user.email={user.email}")
            return redirect(f"/{user.id}")
    else:
        
        form = SignUpForm()
    return render(request, 'registration/user_signup.html', {'form': form})


def user_logout(request):

    logout(request)

    return redirect("/")


@ratelimit(key='ip', rate='10/m', block=True, method=ratelimit.UNSAFE)
def user_login(request):

    if request.method == "POST":
        form = LoginForm(data=request.POST)

        if form.is_valid():
            email = form.cleaned_data['email']
            password = form.cleaned_data['password']

            # Due to an early bug emails may not be unique. Last subscription wins.
            user = User.objects.filter(email__iexact=email).order_by('-id').first()

            if not user:
                form.add_error(None, "This email does not exist.")
                context = dict(form=form)
                return render(request, "registration/user_login.html", context=context)

            user = authenticate(username=user.username, password=password)

            if not user:
                form.add_error(None, "Invalid password.")
            elif user and not user.is_active:
                form.add_error(None, "This user may not log in.")
            elif user and user.is_active:
                login(request, user)
                logger.info(f"logged in user.id={user.id}, user.email={user.email}")

                return redirect(f"/{user.id}")
            else:
                # This should not happen normally.
                form.add_error(None, "Invalid form processing.")
    else:
        initial = dict(nexturl=request.GET.get('next', '/'))
        form = LoginForm(initial)

    context = dict(form=form)
    return render(request, "registration/user_login.html", context=context)


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
    can_create_project = True

    context = dict(projects=projects, steps=steps, can_create_project=can_create_project)

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

    can_create_project = True

    context = dict(projects=project, steps=steps, can_create_project=can_create_project)

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

    can_create_project = True

    context = dict(projects=project, steps=steps, can_create_project=can_create_project,
                   form=form)
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
    can_create_data = True

    context = dict(project=project, steps=steps, can_create_data=can_create_data)

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

    # You can create data from inside another data view
    can_create_data = True

    context = dict(data=data, steps=steps, can_create_data=can_create_data)

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

    project = Project.objects.filter(id=id).first()
    if not project:
        messages.error(request, "No data found for this project.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", False),
        (reverse("analysis_list", kwargs={'id':project.id}),"Analysis List", True),
    ]
    can_create_analysis = True

    context = dict(projects=project, steps=steps, can_create_analysis=can_create_analysis)

    return render(request, "analysis_list.html", context)


#@login_required
def analysis_view(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    if not analysis:
        messages.error(request, f"Data{id} not found.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", False),
        (reverse("analysis_list", kwargs={'id':project.id}),"Analysis List", False),
        (reverse("analysis_view", kwargs={'id': analysis.id}), f"{analysis.title}", True)
    ]

    context = dict(analysis=analysis, steps=steps)

    return render(request, "analysis_view.html", context)


def analysis_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id':project.id}),"Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': analysis.id}), f"{analysis.title}", active)
    ]

    if request.method == "POST":

        form = AnalysisForm(request.POST, instance=analysis)
        if form.is_valid():
            # Redo the dates here to update everytime edited.
            form.save()

    else:
        form = AnalysisForm(instance=analysis)

    context = dict(analysis=analysis, steps=steps, form=form)

    return render(request, 'analysis_edit.html', context)


def analysis_create(request, id):

    project = Project.objects.filter(id=id).first()
    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id': project.id}), "Analysis List", not active),
        (reverse("analysis_create", kwargs={'id': project.id}), "Create analysis", active)
    ]

    if request.method == "POST":

        form = AnalysisForm(data=request.POST)

        if form.is_valid():

            title = form.cleaned_data["title"]
            text = form.cleaned_data["text"]
            owner = User.objects.all().first()
            new_analysis = Analysis.objects.create(title=title, text=text,
                                           owner=owner, project=project)
            new_analysis.save()
            return redirect(reverse("analysis_list", kwargs={'id': project.id}))

        else:
            form.add_error(None, "Invalid form processing.")
    else:

        form = AnalysisForm()
        context = dict(project=project, steps=steps, form=form)
        return render(request, 'analysis_create.html', context)


# also can be seen as results_create
def analysis_run(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project
    # name and label are exepcted in json_file,
    # anything empty will be left
    # Pulls Sequencing Data from needed Project
    #
    json_file = r"""
                [
                {
                    "name": "samples", 
                    "label": "Fastq Samples List", 
                    "selected": "", 
                    "widget": "Select",
                    "form_type": "CharField", 
                    "choices" : {"SE": "single-end", "PE" : "paired-end"}
                },
                {
                    "name": "merge", 
                    "label": "Merge Reads", 
                    "selected": "", 
                    "widget": "Select",
                    "form_type": "CharField",
                    "choices" : {"SE": "single-end", "PE" : "paired-end"}
                },
                ]"""

    active  = True
    owner = User.objects.all().first()

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("analysis_list", kwargs={'id': project.id}), "Analysis List", not active),
        (reverse("analysis_view", kwargs={'id': analysis.id}), f"{analysis.title}", active),
    ]

    if request.method == "POST":
        # loads the json file into the analysis
        #analysis.load(request.POST.json_file)
        form = RunForm(data=request.POST)
        if form.is_valid:

            form.save()
            filled_json = form.json_spec
            filled_makefile = form.makefile

            job = Job.objects.get_or_create(json_data= filled_json,
                                            owner=owner,
                                            analysis=analysis,
                                            makefile=filled_makefile)

            context = dict(job=job, analysis=analysis, steps=steps)
            1/0
            return render(request, "results_list.html", context)

    else:
        # Set json_file in setting.py to avoid loading a file every time a view is served
        analysis.load(json_file)
        print(analysis.json_spec)
        #1/0
        form = RunForm(id=id, json_spec=json_file, makefile=analysis.makefile_template)
        context = dict(analysis=analysis, project=project, steps=steps, form=form)
        return render(request, 'analysis_run.html', context)


def results_list(request, id):

    project = Project.objects.filter(id=id).first()
    analysis_set = project.analysis_set.order_by("-id")

    # Can't populate an empty queryset so populate a list with ids
    # then Results model for those ids
    results = []

    for analysis in analysis_set:

         for result in analysis.result_set.order_by("-id"):
             results.append(result.id)

    # Results belonging to a project
    results = Result.objects.filter(id__in=results)

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", False),
        (reverse("results_list", kwargs={'id': project.id}), "Results List", True),
    ]

    context = dict(project=project, steps=steps, results=results)

    return render(request, "results_list.html", context)


def results_view(request, id):

    result = Result.objects.filter(id=id).first()
    project = result.analysis.project

    if not result:
        messages.error(request, f"Data{id} not found.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", False),
        (reverse("results_list", kwargs={'id':project.id}),"Results List", False),
        (reverse("results_view", kwargs={'id': result.id}), f"{result.title}", True)
    ]

    context = dict(result=result, steps=steps)

    return render(request, "results_view.html", context)


def results_edit(request, id):

    result = Result.objects.filter(id=id).first()
    project = result.analysis.project

    active = True

    steps = [
        (reverse("index"), "Home", not active),
        (reverse("project_list"), "Project List", not active),
        (reverse("project_view", kwargs={'id': project.id}), f"{project.title}", not active),
        (reverse("results_list", kwargs={'id': project.id}), "Results List", not active),
        (reverse("results_view", kwargs={'id': result.id}), f"{result.title}", active)
    ]

    if request.method == "POST":

        form = ResultForm(request.POST, instance=result)
        if form.is_valid():
            # Redo the dates here to update everytime edited.
            form.save()

    else:
        form = ResultForm(instance=result)

    context = dict(result=result, steps=steps, form=form)

    return render(request, 'results_edit.html', context)


def results_create(request, id):

    return











