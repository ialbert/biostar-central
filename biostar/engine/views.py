# import os

import mistune
from django.conf import settings
# from django.template.loader import get_template
from django.utils.safestring import mark_safe
from django.db.models import Q
from django.contrib import messages
from django.views.decorators import cache
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.shortcuts import render, redirect
from django.urls import reverse

from .forms import *
from .const import *
from .decorators import object_access
from .models import (Project, Data, Analysis, Job, User, Access)


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
            step = (reverse("project_view", kwargs={'id': project.id}), PROJECT_ICON, "Project View", is_active)
        elif icon == DATA_LIST_ICON:
            step = (reverse("data_list", kwargs={'id': project.id}), DATA_LIST_ICON, "Data Files", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"File View", is_active)
        elif icon == DATA_UPLOAD:
            step = (reverse("data_view", kwargs={'id': project.id}), DATA_UPLOAD, f"File Upload", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("analysis_list", kwargs={'id': project.id}), ANALYSIS_LIST_ICON, "Analysis Recipes", is_active)
        elif icon == ANALYSIS_VIEW_ICON:
            step = (reverse("analysis_view", kwargs={'id': analysis.id}), ANALYSIS_VIEW_ICON, "Recipe View", is_active)
        elif icon == ANALYSIS_RUN_ICON:
            step = (reverse("analysis_run", kwargs={'id': analysis.id}), ANALYSIS_RUN_ICON, "Analysis Run", is_active)
        elif icon == ANALYSIS_RECIPE_ICON:
            step = (reverse("analysis_recipe", kwargs={'id': analysis.id}), ANALYSIS_RECIPE_ICON, "Recipe Code", is_active)
        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'id': project.id, }), RESULT_LIST_ICON, "Result List", is_active)
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
            step = (reverse("project_view", kwargs={'id': project.id}), ADD_USER, "Manage People", is_active)
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


#TODO: refractor asap
@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def project_users_remove(request, id):
    "Using "
    project = Project.objects.filter(pk=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ADD_USER],
                               project=project)
    searches = []
    access = Access.RECIPE_ACCESS
    # Users with read access to project
    current_users = project.access_set.filter(access__gt=Access.NO_ACCESS)
    current_users = [access.user for access in current_users]

    if request.method == "POST":
        user= request.user

        # Grant users Read access
        form = GrantAccess(data=request.POST, project=project, current_user=user,
                           access=access)
        if form.is_valid(request=request):

            #method = request.POST.get("add_or_remove")
            added, removed, errmsg = form.process(remove=True)
            msg = f"""Removed  <span class="ui green label">{Access.ACCESS_MAP[access]} Permission</span>  
            """
            if errmsg:
                messages.error(request, errmsg)
            else:
                messages.success(request, mark_safe(msg))

        return redirect(reverse("project_users", kwargs=dict(id=project.id)))

    form = GrantAccess(project=project, current_user=request.user, access=access)
    context = dict(steps=steps, current_users=current_users, form=form,
                   available_users=searches, project=project, access=Access(access=access))

    return render(request, "project_users.html", context=context)


#TODO: refractor asap
@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def project_users(request, id):
    "Grants a list of users Read access to a project"
    project = Project.objects.filter(pk=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ADD_USER],
                               project=project)
    searches = []
    access = Access.RECIPE_ACCESS
    # Users with read access to project
    current_users = project.access_set.filter(access__gt=Access.NO_ACCESS)
    current_users = [access.user for access in current_users]

    if request.method == "POST":
        user= request.user
        # Grant users Read access
        form = GrantAccess(data=request.POST, project=project, current_user=user,
                           access=access)
        if form.is_valid(request=request):

            added, removed, errmsg = form.process(add=True)
            msg = f"""Added  <span class="ui green label">{Access.ACCESS_MAP[access]} Permission</span>  
            """
            if errmsg:
                messages.error(request, errmsg)
            else:
                messages.success(request, mark_safe(msg))

        return redirect(reverse("project_users", kwargs=dict(id=project.id)))

    elif request.method == "GET" and request.GET.get("searches"):
        search = request.GET["searches"]
        searches = User.objects.filter( Q(first_name__contains=search) | Q(email__contains=search))
        if not searches:
            messages.info(request, f"No users containing '{search}' found.")

    form = GrantAccess(project=project, current_user=request.user, access=access)
    context = dict(steps=steps, current_users=current_users, form=form,
                   available_users=searches, project=project, access=Access(access=access))

    return render(request, "project_users.html", context=context)


def project_list(request):

    projects = auth.get_project_list(user=request.user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])
    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, id):
    user = request.user
    project = Project.objects.filter(id=id).first()

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

    # Use a placeholder
    access = access or Access(access=Access.NO_ACCESS)
    context = dict(project=project, access=access,
                   data_count=data_count, recipe_count=recipe_count, result_count=result_count,
                   steps=steps)

    return render(request, "project_view.html", context)


@object_access(type=Project, access=Access.EDIT_ACCESS, url='project_view')
def project_edit(request, id):
    project = auth.get_project_list(user=request.user).filter(id=id).first()

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
            owner = request.user

            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream, sticky=sticky, privacy=privacy)
            project.save()

            return redirect(reverse("project_list"))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)
    context = dict(steps=steps, form=form)
    return render(request, 'project_create.html', context)


@object_access(type=Project, access=Access.READ_ACCESS)
def data_list(request, id):

    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON],
                               project=project)
    if not project:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    query = Data.objects.filter(project=project).order_by("sticky", "-date")

    data_list = query.all()
    data_count = query.count()

    context = dict(project=project, steps=steps, data_list=data_list, data_count=data_count)
    return render(request, "data_list.html", context)


# @login_required
@object_access(type=Data, access=Access.READ_ACCESS)
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


@object_access(type=Data, access=Access.EDIT_ACCESS, url='data_view')
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


@object_access(type=Project, access=Access.UPLOAD_ACCESS, url='data_list')
def data_upload(request, id):
    owner = request.user
    project = Project.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_UPLOAD ],
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
            return redirect(reverse("data_list", kwargs={'id': project.id}))
        else:
            form.add_error(None, "Invalid form processing.")
            messages.error(request, "Invalid form processing.")
    else:
        form = DataUploadForm()

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
def analysis_list(request, id):
    """
    Returns the list of analyses for a project id.
    """

    project = Project.objects.filter(id=id).first()
    analyses = Analysis.objects.filter(project=project).order_by("-sticky","-id")

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON],
                               project=project)
    context = dict(project=project, analyses=analyses, steps=steps)

    return render(request, "analysis_list.html", context)


@object_access(type=Analysis, access=Access.READ_ACCESS)
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


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='analysis_view')
def analysis_recipe(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RECIPE_ICON],
                               project=analysis.project, analysis=analysis)

    context = dict(analysis=analysis, steps=steps)
    return render(request, "analysis_recipe.html", context)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='analysis_recipe')
def analysis_copy(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    # Only copy to projects with edit access to
    required_access = Access(access=Access.EDIT_ACCESS)

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RECIPE_ICON],
                               project=analysis.project, analysis=analysis)

    if request.method == "POST":

        if request.user.is_anonymous:
            access_text = Access.ACCESS_MAP.get(Access.RECIPE_ACCESS)
            msg = f"""
            You must be logged in and have the <span class="ui green label">{access_text} Permission</span>  
            to copy an analysis.
            """
            messages.warning(request, mark_safe(msg))
            return redirect(reverse("analysis_copy", kwargs=dict(id=analysis.id)))

        form = AnalysisCopyForm(data=request.POST, analysis=analysis)
        if form.is_valid():
            copied = form.process()
            copied = Project.objects.filter(id=copied[0]).first()
            messages.success(request, f"Copied recipe to {copied.name}.")

        return redirect(reverse("analysis_copy", kwargs=dict(id=analysis.id)))

    projects = auth.get_project_list(user=request.user).exclude(pk=analysis.project.id)
    # Can't touch public projects
    projects = projects.exclude(Q(privacy=Project.PUBLIC))

    # Filter projects by edit access
    if request.user.is_authenticated:
        cond = Q(access__user=request.user, access__access__gt=Access.EDIT_ACCESS)
    else:
        cond = Q(access__access__gt=Access.EDIT_ACCESS)

    projects = projects.filter(cond)

    form = AnalysisCopyForm(analysis=analysis)
    context = dict(analysis=analysis, steps=steps, projects=projects, form=form,
                   project=analysis.project, access=required_access)

    return render(request, "analysis_copy.html", context)

#TODO: refractor asap
@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='analysis_recipe')
def create_copy(request, id):
    "Create a project then copy analysis into it."

    analysis = Analysis.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RECIPE_ICON],
                               project=analysis.project, analysis=analysis)

    if request.user.is_anonymous:
        access_text = Access.ACCESS_MAP.get(Access.RECIPE_ACCESS)
        msg = f"""
        You must be logged in and have the <span class="ui green label">{access_text} Permission</span>  
        to Create Project and Copy Recipe.
        """
        messages.warning(request, mark_safe(msg))
        return redirect(reverse("analysis_copy", kwargs=dict(id=analysis.id)))

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES)

        if form.is_valid():
            name = form.cleaned_data["name"]
            text = form.cleaned_data["text"]
            summary = form.cleaned_data["summary"]
            stream = form.cleaned_data["image"]
            sticky = form.cleaned_data["sticky"]
            privacy = form.cleaned_data["privacy"]
            owner = request.user

            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream, sticky=sticky, privacy=privacy)
            project.save()

            current_params = auth.get_analysis_attr(analysis=analysis, project=project)
            new_analysis = auth.create_analysis(**current_params)
            new_analysis.image.save(analysis.name, analysis.image, save=True)
            new_analysis.name = f"Copy of: {analysis.name}"
            new_analysis.state = analysis.state
            new_analysis.security = analysis.security
            new_analysis.save()
            url = reverse("analysis_list", kwargs=dict(id=new_analysis.project.id))
        else:
            url = reverse("create_copy", kwargs=dict(id=analysis.id))
        return redirect(url)

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)

    context = dict(form=form, steps=steps, analysis=analysis)
    return render(request, "create_copy.html", context)



@object_access(type=Analysis, access=Access.EXECUTE_ACCESS, url='analysis_view')
def analysis_run(request, id):
    analysis = Analysis.objects.filter(id=id).first()

    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
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
    """Function  used to get return updated analysis settings from a given spec"""
    if spec.get("settings"):
        name = spec["settings"].get("name", analysis.name)
        help = spec["settings"].get("help", analysis.text)
        # summary = spec["settings"].get("summary", analysis.text)
        html = make_html(help)

        return dict(name=name, html=html)
    else:
        return dict(name=analysis.name, html=analysis.html)


def process_analysis_edit(analysis, form, method=None):
    form_method_map = {'preview': form.preview,
                       'save': form.save}

    spec = hjson.loads(analysis.json_text)
    json_text = analysis.json_text
    template = analysis.template

    if form.is_valid() and method:

        # Call preview() or save()
        func = form_method_map[method]
        func()
        spec = hjson.loads(form.cleaned_data["json_text"].rstrip())

        # Override json_text and template with most recent
        json_text = form.cleaned_data["json_text"]
        template = form.cleaned_data["template"]

    context = preview_specs(spec, analysis)
    context.update(dict(json_text=json_text,template=template))

    return context


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='analysis_recipe')
def analysis_edit(request, id):
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project
    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON],project=project, analysis=analysis)

    if request.method == "POST":
        form = EditAnalysisForm(analysis=analysis, data=request.POST)
        method = request.POST.get("save_or_preview")
        #Method form.is_valid() called in this function
        context = process_analysis_edit(analysis=analysis, form=form, method=method)
        # should redirect on a save and not a preview

    else:
        form = EditAnalysisForm(analysis=analysis)
        context = process_analysis_edit(analysis=analysis, form=form)

    context.update(dict(project=project, analysis=analysis, steps=steps, form=form))
    return render(request, 'analysis_edit.html', context)


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='analysis_recipe')
def recipe_edit(request, id):

    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

    if request.method == "POST":
        form = AnalysisEditForm(data=request.POST, files=request.FILES, instance=analysis,)
        if form.is_valid():
            form.save()
            return redirect(reverse("analysis_view", kwargs=dict(id=analysis.id)))

    form = AnalysisEditForm(instance=analysis)
    context = dict(steps=steps, analysis=analysis, project=project, form=form)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Project, access=Access.READ_ACCESS)
def job_list(request, id):
    """
    Returns the list of jobs for a project id.
    """
    project = Project.objects.filter(id=id).first()

    if not project:
        messages.error(request, "Jobs not found.")
        # logger.error(f"Jobs for project.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON],
                               project=project)
    jobs =Job.objects.filter(project=project).order_by("-date", "-start_date")

    filter = request.GET.get('filter', '')

    if filter:
        filter = Analysis.objects.filter(id=filter).first()
        jobs = jobs.filter(analysis=filter)

    context = dict(jobs=jobs, steps=steps, project=project, filter=filter)

    return render(request, "job_list.html", context)



@object_access(type=Job, access=Access.EDIT_ACCESS)
def job_edit(request, id):
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON,
                                RESULT_VIEW_ICON],job=job, project=project)

    if request.method == "POST":
        form =JobEditForm(data=request.POST, files=request.FILES, instance=job)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", kwargs=dict(id=job.id)))

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

    context = dict(job=job, steps=steps)
    return render(request, "job_view.html", context=context)


@object_access(type=Job, access=Access.READ_ACCESS)
def job_result_view(request, id):
    """
    Returns the primary result of a job.
    """
    job = Job.objects.filter(id=id).first()
    index = job.json_data.get("settings", {}).get("index", "")

    if job.state == Job.COMPLETED:

        #TODO:This part is still exposed.
        url = settings.MEDIA_URL + job.get_url(path=index)
        return redirect(url)
    else:
        return redirect(reverse("job_view", kwargs=dict(id=id)))


@object_access(type=Job, access=Access.READ_ACCESS)
def job_file_view(request, id):
    """
    Returns the directory view of the job.
    """
    job = Job.objects.filter(id=id).first()
    url = settings.MEDIA_URL + job.get_url()

    return redirect(url)

@object_access(type=Job, access=Access.READ_ACCESS)
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
