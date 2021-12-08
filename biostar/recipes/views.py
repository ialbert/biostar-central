import logging
import os
import toml as hjson
import hashlib
import itertools
import mistune

from django.http import JsonResponse
from django.core.files.storage import default_storage
from django.core.exceptions import PermissionDenied, ImproperlyConfigured
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.db.models import Q, Count
from django.template import loader
from django.db.models import Sum
from django.core.paginator import Paginator
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse
from django.shortcuts import render, redirect, reverse
from django.views.decorators.csrf import ensure_csrf_cookie
from django.template import Template, Context
from django.utils.safestring import mark_safe
from django.core.cache import cache
from ratelimit.decorators import ratelimit
from sendfile import sendfile
from biostar.accounts.models import User
from biostar.recipes import tasks, auth, forms, const, search, util
from biostar.recipes.decorators import read_access, write_access
from biostar.recipes.models import Project, Data, Analysis, Job, Access

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
logger = logging.getLogger('engine')

RATELIMIT_KEY = settings.RATELIMIT_KEY

def join(*args):
    return os.path.abspath(os.path.join(*args))


__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


def index(request):

    # Get the counts from cache
    key = "INDEX_CACHE"

    # Time to live in seconds.
    ttl = 300
    if key not in cache:
        nusers = User.objects.all().count()
        nprojects = Project.objects.all().count()
        nrecipes = Analysis.objects.all().count()
        nresults = Job.objects.all().count()
        value = dict(nusers=nusers, nprojects=nprojects, nrecipes=nrecipes, nresults=nresults)
        cache.set(key, value, ttl)
    else:
        value = cache.get(key, dict())

    context = dict(active="home")
    context.update(value)

    return render(request, 'index.html', context)


@user_passes_test(lambda u: u.is_superuser)
def site_admin(request):
    '''
    Administrative view. Lists the admin project and job.
    '''
    jobs = Job.objects.order_by('-pk')[:200]
    context = dict(jobs=jobs, active="admin")

    return render(request, 'admin_index.html', context=context)


@login_required
def recycle_bin(request):
    "Recycle bin view for a user"
    user = request.user

    BIN_LIMIT = 200

    if user.is_superuser:
        # Super users get access to all deleted objects.
        projects = Project.objects.all()
        query_dict = dict(project__in=projects)
    else:
        # Only searches projects user have access.
        projects = auth.get_project_list(user=user, include_deleted=True)
        query_dict = dict(project__in=projects, owner=user)

    projects = projects.filter(deleted=True).order_by("-lastedit_date")[:BIN_LIMIT]

    # Filter data, recipes, and jobs according to projects user has access to.
    data = Data.objects.filter(**query_dict, deleted=True).order_by("-lastedit_date")[:BIN_LIMIT]

    recipes = Analysis.objects.filter(**query_dict, deleted=True).order_by("-lastedit_date")[:BIN_LIMIT]

    jobs = Job.objects.filter(**query_dict, deleted=True).order_by("-lastedit_date")[:BIN_LIMIT]

    deleted = []

    for obj in [projects, data, recipes, jobs]:
        for item in obj:
            deleted.append(item)

    deleted = sorted(deleted, key=lambda x: x.lastedit_date, reverse=True)
    context = dict(deleted=deleted, active="bin")

    return render(request, 'recycle_bin.html', context=context)


@write_access(type=Project, fallback_view="project_view")
def project_delete(request, uid):
    project = Project.objects.filter(uid=uid).first()
    project.deleted = not project.deleted
    project.save()

    # Same function called tor restore projects.
    if project.deleted:
        msg = f"Project:{project.name} deleted"
    else:
        msg = f"Project:{project.name} restored"

    messages.success(request, msg)

    return redirect(reverse("project_list"))


def search_bar(request):
    results = search.search(request=request)

    # Indicate to users that minimum character needs to be met.
    query_lenth = len(request.GET.get("q", "").strip())
    min_length = query_lenth > settings.SEARCH_CHAR_MIN

    # Indicate to users that there are no results for search.
    no_results = min_length and len(results) == 0

    context = dict(results=results, query=request.GET.get("q", "").strip(),
                   min_length=min_length, no_results=no_results)

    return render(request, "search.html", context)


@write_access(type=Project, fallback_view="data_list")
def project_users(request, uid):
    """
    Manage project users page
    """
    project = Project.objects.filter(uid=uid).first()
    # Get users that already have access to project.
    have_access = project.access_set.exclude(access=Access.NO_ACCESS).order_by('-date')

    # Search query for users.
    q = request.GET.get("q", "")
    if q:
        selected = have_access.values('user')
        targets = User.objects.filter(Q(email__contains=q) | Q(profile__name__contains=q) |
                                      Q(username__contains=q) |
                                      Q(profile__uid__contains=q)).exclude(id__in=selected)
    else:
        targets = []

    # Gather access forms for users who currently have access
    context = dict(project=project, have_access=have_access, q=q, targets=targets, activate='User Management')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "project_users.html", context=context)


@read_access(type=Project)
def project_info(request, uid):
    user = request.user

    project = Project.objects.filter(uid=uid).first()

    # Show counts for the project.
    counts = get_counts(project)

    # Who has write access
    write_access = auth.is_writable(user=user, project=project)
    if user.is_authenticated:
        access = Access.objects.filter(user=user, project=project).first()
    else:
        access = Access(access=Access.NO_ACCESS)

    access = access or Access(access=Access.NO_ACCESS)

    context = dict(project=project, active="info", write_access=write_access, access=access)
    context.update(counts)

    return render(request, "project_info.html", context)


def project_list(request):
    user = request.user
    projects = auth.get_project_list(user=user)
    page = request.GET.get("page")
    projects = projects.order_by("-rank")

    # Add pagination.
    paginator = Paginator(projects, per_page=settings.PER_PAGE)
    projects = paginator.get_page(page)

    context = dict(projects=projects, active="project_list")
    return render(request, "project_list.html", context=context)


def latest_recipes(request):
    """
    """
    page = request.GET.get("page")
    # Select public recipes
    recipes = Analysis.objects.filter(project__privacy=Project.PUBLIC, deleted=False)
    recipes = recipes.order_by("-rank", "-lastedit_date")[:50]

    recipes = recipes.annotate(job_count=Count("job", filter=Q(job__deleted=False)))

    paginator = Paginator(recipes, per_page=settings.PER_PAGE)
    recipes = paginator.get_page(page)

    context = dict(recipes=recipes, active="latest_recipes")

    return render(request, "latest_recipes.html", context=context)


@read_access(type=Project)
def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="data_list.html",
                        active='data', show_summary=True)


@read_access(type=Project)
def recipe_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes')


def job_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


def get_counts(project, user=None):
    data_count = project.data_count
    recipe_count = project.recipes_count

    result_count = project.jobs_count
    discussion_count = 0

    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count,
        discussion_count=discussion_count
    )


@read_access(type=Project)
def project_view(request, uid, template_name="project_info.html", active='info', show_summary=None,
                 extra_context={}):
    """
    This view handles the project info, data list, recipe list, result list views.
    """
    page = request.GET.get('page')

    # The user making the request
    user = request.user

    # The project that is viewed.
    project = Project.objects.filter(uid=uid).first()

    # Select all the data in the project.
    data_list = project.data_set.filter(deleted=False).order_by("-lastedit_date", "rank", "-date").all()
    data_paginator = Paginator(data_list, per_page=settings.PER_PAGE)
    data_list = data_paginator.get_page(page)

    recipe_list = project.analysis_set.filter(deleted=False).order_by("-rank", "-lastedit_date", "-date").all()

    # Annotate each recipe with the number of jobs it has.
    recipe_list = recipe_list.annotate(job_count=Count("job", filter=Q(job__deleted=False)))
    recipe_paginator = Paginator(recipe_list, per_page=settings.PER_PAGE)
    recipe_list = recipe_paginator.get_page(page)

    job_list = project.job_set.filter(deleted=False).order_by("-lastedit_date").all()

    # Filter job results by analysis
    filter_uid = request.GET.get('filter', '')
    recipe_filter = Analysis.objects.filter(uid=filter_uid).first()

    # The recipe filter exists
    if recipe_filter:
        job_list = job_list.filter(analysis=recipe_filter)

    # Add related content.
    job_list = job_list.select_related("analysis")
    job_paginator = Paginator(job_list, per_page=settings.PER_PAGE)
    job_list = job_paginator.get_page(page)

    # Who has write access
    write_access = auth.is_writable(user=user, project=project)

    # Build the context for the project.
    context = dict(project=project, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, recipe_filter=recipe_filter, write_access=write_access, rerun_btn=True,
                   include_copy=False)

    # Compute counts for the project.
    counts = get_counts(project)

    # Update conext with the counts.
    context.update(counts)

    # Add any extra context that may come from parameters.
    context.update(extra_context)

    return render(request, template_name, context)


@write_access(type=Project, fallback_view="data_list")
def project_edit(request, uid):
    "Edit meta-data associated with a project."

    project = Project.objects.filter(uid=uid).first()
    form = forms.ProjectForm(instance=project, request=request)
    if request.method == "POST":
        form = forms.ProjectForm(data=request.POST, files=request.FILES, instance=project, request=request)
        if form.is_valid():
            project = form.save()
            Project.objects.filter(uid=uid).update(lastedit_user=request.user)
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

    context = dict(project=project, form=form, activate='Edit Project')

    context.update(get_counts(project))
    return render(request, "project_edit.html", context=context)


@login_required
@ratelimit(key=RATELIMIT_KEY, rate='5/h', block=True, method=ratelimit.UNSAFE)
def project_create(request):
    """
    View used create an empty project belonging to request.user.
    Input is validated with a form and actual creation is routed through auth.create_project.
    """
    user = request.user
    project = auth.create_project(user=user, text="Project information goes here. ", name="Project")

    # Ensure project name is unique.
    project.name = f"{project.name} {project.id}"
    project.text = f"Add more information on {project.name}"
    project.save()

    messages.success(request, "A new project has been created. Please set the project details")
    return redirect(reverse("project_edit", kwargs=dict(uid=project.uid)))


@read_access(type=Data)
def data_view(request, uid):
    "Show information specific to each data."

    data = Data.objects.filter(uid=uid).first()
    project = data.project
    paths = auth.listing(root=data.get_data_dir())

    context = dict(data=data, project=project, paths=paths, serve_view="data_serve",
                   activate='Selected Data', uid=data.uid, show_all=True)
    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@write_access(type=Data, fallback_view="data_view")
def data_edit(request, uid):
    """
    Edit meta-data associated with Data.
    """
    data = Data.objects.filter(uid=uid).first()
    form = forms.DataEditForm(instance=data, initial=dict(type=data.type), user=request.user)

    if request.method == "POST":
        form = forms.DataEditForm(data=request.POST, instance=data, user=request.user, files=request.FILES)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(uid=data.uid)))

    context = dict(data=data, form=form, activate='Edit Data', project=data.project)

    context.update(get_counts(data.project))
    return render(request, 'data_edit.html', context)


@write_access(type=Project, fallback_view="data_list")
def data_upload(request, uid):
    "Data upload view routed through auth.create_data."
    logger.info("Entering data upload view.")

    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    form = forms.DataUploadForm(user=owner, project=project)

    if request.method == "POST":
        logger.info("Initating data upload with POST.")
        form = forms.DataUploadForm(data=request.POST, files=request.FILES, user=owner, project=project)
        if form.is_valid():
            data = form.save()
            logger.info("Data upload complete, redirecting.")
            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))
        
    logger.info(f"Computing context values.")
    uploaded_files = Data.objects.filter(owner=owner, method=Data.UPLOAD)
    # The current size of the existing data
    current_size = uploaded_files.aggregate(Sum("size"))["size__sum"] or 0
    # Maximum data that may be uploaded.
    maximum_size = owner.profile.upload_size * 1024 * 1024

    context = dict(project=project, form=form,
                   maximum_size=maximum_size, activate='Upload data',
                   current_size=current_size)

    counts = get_counts(project)
    context.update(counts)
    
    logger.info(f"Rendering data upload template. ")
    return render(request, 'data_upload.html', context)


@read_access(type=Analysis)
def recipe_code_download(request, uid, fname):
    """
    Download the raw recipe template as a file
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    script = auth.render_script(recipe=recipe)

    # Trigger file download with given filename.
    response = HttpResponse(script, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename={fname}'

    return response


@read_access(type=Analysis)
@ratelimit(key=RATELIMIT_KEY, rate='10/h', block=True, method=ratelimit.UNSAFE)
def recipe_run(request, uid):
    """
    View used to execute recipes and start a 'Queued' job.
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    # Form submission.
    if request.method == "POST":
        form = forms.RecipeInterface(request=request, analysis=recipe, json_data=recipe.json_data,
                                     data=request.POST, files=request.FILES)
        # The form validation will authorize the job.
        if form.is_valid():
            # Create the job from the recipe and incoming json data.
            job = auth.create_job(analysis=recipe, user=request.user, fill_with=form.cleaned_data)
            # Spool via UWSGI or start it synchronously.
            tasks.execute_job.spool(job_id=job.id)
            url = reverse("recipe_view", kwargs=dict(uid=recipe.uid)) + "#results"
            return redirect(url)
        else:
            messages.error(request, form.errors)

    url = reverse("recipe_view", kwargs=dict(uid=recipe.uid)) + "#run"

    return redirect(url)


@read_access(type=Job)
def job_rerun(request, uid):
    # Get the job.
    job = Job.objects.filter(uid=uid).first()
    next = request.GET.get('next')
    # Get the recipe
    recipe = job.analysis
    # Get the job JSON
    json_data = job.json_data

    # Validate users can run the recipe.
    valid, msg = auth.validate_recipe_run(user=request.user, recipe=recipe)

    if not valid:
        messages.error(request, msg)
        redir = next or reverse('job_view', kwargs=dict(uid=job.uid))
        return redirect(redir)

    # Create a new job
    job = auth.create_job(analysis=recipe, user=request.user, json_data=json_data)

    # Spool via UWSGI or run it synchronously.
    tasks.execute_job.spool(job_id=job.id)
    if auth.is_readable(user=request.user, obj=recipe):
        url = reverse('recipe_view', kwargs=dict(uid=job.analysis.uid)) + "#results"
    else:
        url = job.url()

    return redirect(url)


def get_part(request, name, id):
    """
    Return a template by name and with uid rendering
    """

    user = request.user

    # The recipe that needs to be edited.
    recipe = Analysis.objects.filter(id=id).annotate(
        job_count=Count("job", filter=Q(job__deleted=False))
    ).first()

    project = recipe.project

    if not auth.is_readable(obj=project, user=user):
        message = str("Recipe is not readable by current user")
        return HttpResponse(message)

    # Fills in project level counts (results, data and recipe counts).
    counts = get_counts(recipe.project)

    if name == "run":
        initial = dict(name=f"Results for: {recipe.name}")
        form = forms.RecipeInterface(request=request, analysis=recipe,
                                     json_data=recipe.json_data, initial=initial)
    else:
        # Initial form loading via a GET request.
        form = forms.RecipeForm(instance=recipe, user=request.user, project=project)


    remap = dict(
        info="parts/recipe_info.html",
        code="parts/recipe_code.html",
        interface="parts/recipe_interface.html",
        run="parts/recipe_run.html",
        results="parts/recipe_results.html",
        details='parts/recipe_details.html',
    )

    name = remap.get(name, "parts/placeholder.html")

    # Check to see if this recipe is runnable by the user.
    is_runnable = auth.authorize_run(user=user, recipe=recipe)

    # Check to see if recipe is editable
    editable = auth.writeable_recipe(user=user, source=recipe)

    # Get the list of jobs required for recipe results
    jobs = recipe.job_set.filter(deleted=False).order_by("-lastedit_date").all()
    context = dict(recipe=recipe, form=form, is_runnable=is_runnable, job_list=jobs, rerun_btn=False,
                   include_copy=False, editable=editable, user=user, project=recipe.project)
    context.update(counts)

    html = render(request, name, context=context)
    return html


@ensure_csrf_cookie
@read_access(type=Analysis)
def recipe_view(request, uid):
    """
    Edit meta-data associated with a recipe.
    """

    # The user making the request.
    user = request.user

    # The recipe that needs to be edited.
    recipe = Analysis.objects.filter(uid=uid).annotate(
            job_count=Count("job", filter=Q(job__deleted=False))
        ).first()

    # The project that recipe belongs to.
    project = recipe.project

    # Initial form loading via a GET request.
    form = forms.RecipeForm(instance=recipe, user=request.user, project=project)

    # Fills in project level counts (results, data and recipe counts).
    counts = get_counts(project)

    # Disable buttons if project not writeable.
    btn_state = '' if auth.is_writable(user=user, project=project) else 'disabled'

    # Get the list of jobs required to view recipe results
    jobs = recipe.job_set.filter(deleted=False).order_by("-lastedit_date").all()

    # Check to see if this recipe is runnable by the user.
    is_runnable = auth.authorize_run(user=user, recipe=recipe)

    # Check to see if recipe is editable
    editable = auth.writeable_recipe(user=user, source=recipe)

    # Generate the context.
    context = dict(recipe=recipe, job_list=jobs, project=project, form=form, btn_state=btn_state,
                   is_runnable=is_runnable, activate='Recipe View', rerun_btn=False,
                   include_copy=False, editable=editable)

    # Update context with counts.
    context.update(counts)

    return render(request, 'recipe_view.html', context)


@read_access(type=Project, strict=True)
def recipe_create(request, uid):
    # Get the project
    project = Project.objects.filter(uid=uid).first()
    recipe = auth.create_analysis(project=project, name="Recipe", template="echo 'Hello World!'")

    # Ensure recipe names are distinguishable from one another.
    recipe.name = f"{recipe.name} {recipe.id}"
    recipe.save()

    url = reverse("recipe_view", kwargs=dict(uid=recipe.uid))
    messages.success(request, "A new recipe has been created.")
    return redirect(url)


@write_access(type=Job, fallback_view="job_view")
def job_edit(request, uid):
    "Edit meta-data associated with a job."

    job = Job.objects.filter(uid=uid).first()
    project = job.project
    form = forms.JobEditForm(instance=job, user=request.user)

    if request.method == "POST":
        form = forms.JobEditForm(data=request.POST, files=request.FILES, instance=job, user=request.user)
        if form.is_valid():
            form.save()
            return redirect(reverse("job_view", kwargs=dict(uid=job.uid)))

    context = dict(job=job, project=project, form=form, activate='Edit Result')

    context.update(get_counts(project))
    return render(request, 'job_edit.html', context)


@write_access(type=Analysis, fallback_view="recipe_view")
def recipe_delete(request, uid):
    recipe = Analysis.objects.filter(uid=uid).first()
    user = request.user

    auth.delete_recipe(recipe=recipe, user=user)
    tmpl = loader.get_template('widgets/delete_msg.html')
    context = dict(obj=recipe, undo_url=reverse('recipe_delete', kwargs=dict(uid=recipe.uid)))
    msg = tmpl.render(context=context)

    messages.success(request, mark_safe(msg))

    return redirect(reverse("recipe_list", kwargs=dict(uid=recipe.project.uid)))


@write_access(type=Job, fallback_view="job_view")
def job_delete(request, uid):
    job = Job.objects.filter(uid=uid).first()

    running_job = job.state == Job.RUNNING and not job.deleted

    if running_job:
        messages.error(request, "Can not delete a running job. Wait until it finishes.")
        return redirect(job.url())

    auth.delete_object(obj=job, request=request)
    msg = f"Deleted <b>{job.name}</b>." if job.deleted else f"Restored <b>{job.name}</b>."
    messages.success(request, mark_safe(msg))

    return redirect(reverse("job_list", kwargs=dict(uid=job.project.uid)))


@write_access(type=Data, fallback_view="data_view")
def data_delete(request, uid):
    data = Data.objects.filter(uid=uid).first()

    auth.delete_object(obj=data, request=request)
    msg = f"Deleted <b>{data.name}</b>." if data.deleted else f"Restored <b>{data.name}</b>."
    messages.success(request, mark_safe(msg))

    return redirect(reverse("data_list", kwargs=dict(uid=data.project.uid)))


@read_access(type=Job)
def job_view(request, uid):
    '''
    Views the state of a single job.
    '''
    job = Job.objects.filter(uid=uid).first()
    project = job.project

    stdout = job.stdout_log
    stderr = job.stderr_log
    if job.is_running():
        # Pass the most current stderr and stdout
        stdout_path = os.path.join(job.path, settings.JOB_STDOUT)
        stderr_path = os.path.join(job.path, settings.JOB_STDERR)
        stdout = open(stdout_path, 'r').read()
        stderr = open(stderr_path, 'r').read()

    paths = auth.listing(root=job.get_data_dir())

    # Pass along any plugins this job has.
    plugin = job.json_data.get('settings', {}).get('plugin')

    context = dict(job=job, project=project, stderr=stderr, stdout=stdout,uid=job.uid, show_all=True,
                   activate='View Result', paths=paths, serve_view="job_serve",
                   plugin=plugin)

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def file_serve(request, path, obj, download=False):
    """
    Authenticates access through decorator before serving file.
    """

    # Get the object that corresponds to the entry.
    root = obj.get_data_dir()

    # This will turn into an absolute path.
    file_path = join(root, path)

    # Ensure only files in the object root can be accessed.
    if not file_path.startswith(root) or not os.path.isfile(file_path):
        msg = "Invalid path." if (not file_path.startswith(root)) else f"File not found: {path}"
        messages.error(request, msg)
        return redirect(obj.url())

    # The response will be the file content.
    mimetype = auth.guess_mimetype(fname=path)

    # Get the filesize in Mb
    size = os.path.getsize(file_path) / 1024 / 1024

    # This behavior can be further customized in front end webserver.
    if size < 20 and not download:
        # Return small files in the browser if possible.
        data = sendfile(request, file_path, mimetype=mimetype)
    else:
        # Trigger a file download for bigger files.
        fname = os.path.basename(file_path)
        data = sendfile(request, file_path, attachment=True, attachment_filename=fname, mimetype=mimetype)

    return data


@read_access(type=Data, allowed_cors='view.qiime2.org')
def data_serve(request, uid, path):
    """
    Serves files from a data directory.
    """
    obj = Data.objects.filter(uid=uid).first()
    return file_serve(request=request, path=path, obj=obj)


@read_access(type=Data)
def data_download(request, uid):
    """
    Download a given data object.
    """
    obj = Data.objects.filter(uid=uid).first()
    files = obj.get_files()
    # Get first file in data directory.
    path = files[0]

    return file_serve(request=request, path=path, obj=obj, download=True)


def job_serve(request, uid, path):
    """
    Serves files from a job directory.
    """
    obj = Job.objects.filter(uid=uid).first()

    if obj:
        return file_serve(request=request, path=path, obj=obj)
    else:
        messages.error(request, "Object does not exist")
        return redirect("/")


@login_required
def project_share(request, token):
    # Get the project by the it's unique share token.
    project = Project.objects.filter(sharable_token=token).first()
    user = request.user

    if not project.is_shareable:
        messages.error(request, "Project is not sharable.")
        return redirect(reverse('project_list'))

    access = Access.objects.filter(user=user, project=project).first()
    # Give user Share access.xs

    # Create access
    if access is None:
        Access.objects.create(user=user, project=project, access=Access.SHARE_ACCESS)

    # Update existing No Access to Share
    elif access.access == Access.NO_ACCESS:
        Access.objects.filter(id=access.id).update(access=Access.SHARE_ACCESS)

    messages.success(request, "Granted share access")
    return redirect(reverse('project_view', kwargs=dict(uid=project.uid)))


@login_required
def import_files(request, path=""):
    """
    Import files mounted on IMPORT_ROOT_DIR in settings
    """
    user = request.user
    if not user.profile.trusted:
        messages.error(request, 'Only trusted users may views this page.')
        return redirect(reverse('project_list'))

    root = settings.IMPORT_ROOT_DIR
    path = os.path.abspath(os.path.join(root, path))
    if not path.startswith(root):
        messages.error(request, 'Outside root directory')
        path = ''

    if not os.path.exists(path):
        messages.error(request, 'File path does not exist')
        path = ''

    # Set the current node we are traversing.
    node = path if os.path.isdir(path) else None

    # Walk through the /root/node/ and collect paths.
    # Directories are not walked through because show_all=False.
    paths = auth.listing(root=root, node=node, show_all=False)

    context = dict(paths=paths, active="import", show_all=False)

    return render(request, 'import_files.html', context=context)


