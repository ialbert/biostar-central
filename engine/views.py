import uuid
import logging

from django.contrib.auth.decorators import login_required
from django.contrib.auth import login, authenticate, logout
from django.contrib import messages
from django.shortcuts import render, redirect
from django.urls import reverse

from .forms import SignUpForm, LoginForm, ProjectForm, DataForm

from .models import (User, Project, Data,
                     Analysis)
from ratelimit.decorators import ratelimit

logger = logging.getLogger('engine')


def index(request):

    steps = [
        (reverse("index"), "Home", True)
    ]
    context = dict(steps=steps)
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


@login_required
def home(request, id):
    return render(request, 'home.html')


#@login_required
def project_list(request):

    projects = Project.objects.all()
    if not projects:
        messages.error(request, "No project found.")
        return redirect("/")

    # True for the active section at the moment
    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", True)
    ]
    is_project_list = True

    context = dict(projects=projects, steps=steps, is_project_list=is_project_list)

    return render(request, "project_list.html", context)


#@login_required
def project_view(request, id):

    project = Project.objects.filter(id=id).first()
    if not project:
        messages.error(request, f"Project{id} not found.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view", kwargs={'id':project.id}), f"{project.title}", True)
    ]
    is_project_list = True

    context = dict(projects=project, steps=steps, is_project_list=is_project_list)

    return render(request, "project_view.html", context)


def project_create(request):

    #else:
    #form = ProjectForm()
    return


def project_edit(request, id):

    if request.method == "POST":

        proj = Project.objects.filter(id=id).first()
        form = ProjectForm(request.POST, instance=proj)

        if form.is_valid():
            form.save()

            return render(request, 'project/project_edit.html',
                          {'form': form, 'object_list': proj})

    else:
        proj = Project.objects.filter(id=id).first()
        form = ProjectForm(instance=proj)

        return render(request, 'project/project_edit.html',
                      {'form': form, 'object_list': proj})


#@login_required
def data_list(request, id):

    project = Project.objects.filter(id=id).first()
    if not project:
        messages.error(request, "No data found for this project.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", False),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", True),
    ]
    is_data_list = True

    context = dict(projects=project, steps=steps, is_data_list=is_data_list)

    return render(request, "data_list.html", context)


#@login_required
def data_view(request, id):

    data = Data.objects.filter(id=id).first()
    project = data.project

    if not data:
        messages.error(request, f"Data{id} not found.")

    steps = [
        (reverse("index"), "Home", False),
        (reverse("project_list"), "Project List", False),
        (reverse("project_view",kwargs={'id':project.id}), f"{project.title}", False),
        (reverse("data_list", kwargs={'id':project.id}),"Data List", False),
        (reverse("data_view", kwargs={'id': data.id}), f"{data.title}", True)
    ]

    # You can create data from inside another data view
    can_create = True

    context = dict(data=data, steps=steps, is_data_list=can_create)

    return render(request, "data_view.html", context)


def data_create(request, id):

    proj = Project.objects.filter(id=id).first()
    created_data_id = Data.objects.filter(project=proj).values_list('id', flat=True)

    if request.method == "POST":

        current_data = created_data_id[len(created_data_id)-1]
        data = Data.objects.filter(id=current_data).first()

        form = DataForm(data=request.POST, instance=data)

        if form.is_valid():
            form.save()

            return render(request, 'project/data_create.html',
                          {'form': form, 'object_list': proj})
        else:
            print(form.is_valid)
            print(form.cleaned_data["title"])
            print(form.cleaned_data["owner"])


            return render(request, 'project/data_create.html',
                          {'form': form, 'object_list': proj})
    else:

        form = DataForm()

        return render(request, 'project/data_create.html',
                      {'form': form, 'object_list': proj})


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
    is_analysis_list = True

    context = dict(projects=project, steps=steps, is_analysis_list=is_analysis_list)

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

    can_create = True

    context = dict(analysis=analysis, steps=steps, is_analysis_list=can_create)

    return render(request, "analysis_view.html", context)

def analysis_list_edit(request, id):

    return None


def analysis_detail_edit(request, id):

    return None
