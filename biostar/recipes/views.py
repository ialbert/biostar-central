import logging
import os
import toml as hjson
import hashlib
import itertools
import mistune
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.decorators import user_passes_test
from django.db.models import Q, Count
from django.db.models import Sum
from django.http import HttpResponse
from django.shortcuts import render, redirect, reverse
from django.template import Template, Context
from django.utils.safestring import mark_safe
from ratelimit.decorators import ratelimit
from sendfile import sendfile
from biostar.accounts.models import User
from biostar.recipes import tasks, auth, forms, const, search, util
from biostar.recipes.decorators import read_access, write_access, exists
from biostar.recipes.models import Project, Data, Analysis, Job, Access

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
logger = logging.getLogger('engine')


def join(*args):
    return os.path.abspath(os.path.join(*args))


__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


def index(request):
    context = dict(active="home")
    return render(request, 'index.html', context)



@user_passes_test(lambda u: u.is_superuser)
def site_admin(request):
    '''
    Administrative view. Lists the admin project and job.
    '''
    jobs = Job.objects.order_by('-pk')[:200]
    context = dict(jobs=jobs)

    return render(request, 'admin_index.html', context=context)


def about(request):
    """
    Added an about page with the
    """

    # Get the docs
    try:
        recipe_docs = os.path.join(settings.DOCS_ROOT, 'recipes', 'recipes.md')
        recipe_docs = open(recipe_docs, 'r').read()
        html = mistune.markdown(recipe_docs, escape=False)
    except Exception as exc:
        logger.error(f'Error loading about page: {exc}')
        html = "About page"
    html = mark_safe(html)
    context = dict(html=html)
    return render(request, 'about.html', context=context)


@login_required
def recycle_bin(request):
    "Recycle bin view for a user"
    user = request.user

    if user.is_superuser:
        # Super users get access to all deleted objects.
        projects = Project.objects.all()
        query_dict = dict(project__in=projects)
    else:
        # Only searches projects user have access.
        projects = auth.get_project_list(user=user, include_deleted=True)
        query_dict = dict(project__in=projects, owner=user)

    projects = projects.filter(deleted=True).order_by("date")
    data = Data.objects.filter(**query_dict, deleted=True).order_by("date")
    recipes = Analysis.objects.filter(**query_dict, deleted=True).order_by("date")
    jobs = Job.objects.filter(**query_dict, deleted=True).order_by("date")

    context = dict(jobs=jobs, data=data, recipes=recipes, projects=projects, active="bin")

    return render(request, 'recycle_bin.html', context=context)


@write_access(type=Project, fallback_view="project_view")
def project_delete(request, uid):
    project = Project.objects.filter(uid=uid).first()
    project.deleted = not project.deleted
    project.save()

    msg = f"Project:{project.name} successfully "
    msg += "deleted!" if project.deleted else "restored!"

    messages.success(request, msg)

    return redirect(reverse("project_list_private"))


def clear_clipboard(request, uid):
    "Clear copied objects held in clipboard."

    next_url = request.GET.get("next", reverse("project_view", kwargs=dict(uid=uid)))
    board = request.GET.get("board")
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})

    if clipboard.get(board):
        clipboard[board] = []
        request.session.update({settings.CLIPBOARD_NAME: clipboard})

    return redirect(next_url)


def search_bar(request):
    results = search.search(request=request)
    # Indicate to users that minimum character needs to be met.
    query_lenth = len(request.GET.get("q", "").strip())
    min_length = query_lenth > settings.SEARCH_CHAR_MIN

    # Indicate to users that there are no results for search.
    current_results = len([inner for outer in results.values() for inner in outer])
    no_results = min_length and current_results == 0

    context = dict(results=results, query=request.GET.get("q", "").strip(),
                   min_length=min_length, no_results=no_results)

    return render(request, "search.html", context)


def data_download(request, uid):
    """
    Download data of given uid.
    """

    return


@write_access(type=Project, fallback_view="data_list")
def project_users(request, uid):
    """
    Manage project users page
    """
    project = Project.objects.filter(label=uid).first()
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


def project_list_private(request):
    """Only list private projects belonging to a user."""

    projects = auth.get_project_list(user=request.user, include_public=False)

    empty_msg = "No projects found."
    if request.user.is_anonymous:
        projects = []
        empty_msg = mark_safe(f"You need to <a href={reverse('login')}> log in</a> to view your projects.")
    else:
        projects = projects.order_by("rank", "-date", "-lastedit_date", "-id")

    context = dict(projects=projects, empty_msg=empty_msg, active="projects", icon='briefcase',
                   title='Private Projects',
                   private='active')

    return render(request, "project_list.html", context)


def project_list_public(request):
    """Only list public projects."""

    projects = auth.get_project_list(user=request.user)
    # Exclude private projects
    projects = projects.exclude(privacy__in=[Project.PRIVATE, Project.SHAREABLE])
    projects = projects.order_by("rank", "-date", "-lastedit_date", "-id")

    context = dict(projects=projects, active="projects", icon='list', title='Public Projects',
                   public='active', empty_msg="No projects found.")

    return render(request, "project_list.html", context)


def project_list(request):
    if request.user.is_authenticated:
        # Return private projects when user is logged in.
        return project_list_private(request)
    else:
        return project_list_public(request)


@read_access(type=Project)
def data_list(request, uid):
    """
    Returns the list of data for a project uid.
    """
    extra_context = dict(copied_data=const.COPIED_DATA,
                         data_paste_targets=const.DATA_PASTE_TARGETS)
    return project_view(request=request, uid=uid, template_name="data_list.html",
                        active='data', show_summary=True, extra_context=extra_context)


@read_access(type=Project)
def recipe_list(request, uid):
    """
    Returns the list of recipes for a project uid.
    """
    extra_context = dict(recipe_paste_targets=const.COPIED_RECIPES)

    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes',
                        extra_context=extra_context)


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


@exists(otype=Project)
def recipe_listing(request, label):
    project = Project.objects.filter(label=label).first()
    return recipe_list(request, uid=project.uid)


@exists(otype=Project)
def project_viewing(request, label):
    project = Project.objects.filter(label=label).first()
    print(label, "FOOOOOOOO")
    return project_view(request=request, uid=project.uid)


@exists(otype=Project)
def job_listing(request, label):
    project = Project.objects.filter(label=label).first()
    return job_list(request=request, uid=project.uid)


@exists(otype=Project)
def data_listing(request, label):
    project = Project.objects.filter(label=label).first()
    return data_list(request, uid=project.uid)


@read_access(type=Project)
def project_view(request, uid, template_name="project_info.html", active='info', show_summary=None,
                 extra_context={}):
    """
    This view handles the project info, data list, recipe list, result list views.
    """

    # The user making the request
    user = request.user

    # The project that is viewed.
    project = Project.objects.filter(uid=uid).first()

    # Select all the data in the project.
    data_list = project.data_set.filter(deleted=False).order_by("rank", "-date").all()
    recipe_list = project.analysis_set.filter(deleted=False).order_by("rank", "-date").all()

    # Annotate each recipe with the number of jobs it has.
    recipe_list = recipe_list.annotate(job_count=Count("job", filter=Q(job__deleted=False)))

    job_list = project.job_set.filter(deleted=False).order_by("-lastedit_date").all()

    # Filter job results by analysis
    filter_uid = request.GET.get('filter', '')
    recipe_filter = Analysis.objects.filter(uid=filter_uid).first()

    # The recipe filter exists
    if recipe_filter:
        job_list = job_list.filter(analysis=recipe_filter)

    # Add related content.
    job_list = job_list.select_related("analysis")

    # Who has write access
    write_access = auth.is_writable(user=user, project=project)

    # Build the context for the project.
    context = dict(project=project, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, recipe_filter=recipe_filter, write_access=write_access)

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
@ratelimit(key='ip', rate='5/h', block=True, method=ratelimit.UNSAFE)
def project_create(request):
    """
    View used create an empty project belonging to request.user.
    Input is validated with a form and actual creation is routed through auth.create_project.
    """
    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = forms.ProjectForm(initial=initial, request=request, create=True)

    if request.method == "POST":
        # create new projects here ( just populates metadata ).
        form = forms.ProjectForm(request=request, data=request.POST, create=True, files=request.FILES)
        if form.is_valid():
            project = form.custom_save(owner=request.user)
            print(project, "LOLOLLLL", project.label)

            return redirect(reverse("project_viewing", kwargs=dict(label=project.label)))

    context = dict(form=form)
    return render(request, "project_create.html", context=context)


@read_access(type=Data)
def data_copy(request, uid):
    data = Data.objects.filter(uid=uid).first()
    next_url = request.GET.get("next", reverse("data_list", kwargs=dict(uid=data.project.uid)))

    board_items = auth.copy_uid(request=request, uid=data.uid, board=const.COPIED_DATA)
    messages.success(request, f"Copied item(s), clipboard contains {len(set(board_items))}.")
    return redirect(next_url)


@read_access(type=Analysis)
def recipe_copy(request, uid):
    recipe = Analysis.objects.filter(uid=uid).first()
    next_url = request.GET.get("next", reverse("recipe_list", kwargs=dict(uid=recipe.project.uid)))

    board_items = auth.copy_uid(request=request, uid=recipe.uid, board=const.COPIED_RECIPES)
    messages.success(request, f"Copied recipe, you currently have  {len(set(board_items))} copied.")
    return redirect(next_url)


@read_access(type=Job)
def job_copy(request, uid):
    job = Job.objects.filter(uid=uid).first()
    next_url = request.GET.get("next", reverse("job_list", kwargs=dict(uid=job.project.uid)))

    board_items = auth.copy_uid(request=request, uid=job.uid, board=const.COPIED_RESULTS)
    messages.success(request, f"Copied item(s), clipboard contains {len(set(board_items))}.")
    return redirect(next_url)


@read_access(type=Data)
def data_file_copy(request, uid, path):
    # Get the root data where the file exists
    data = Data.objects.filter(uid=uid).first()
    fullpath = os.path.join(data.get_data_dir(), path)
    copied = auth.copy_file(request=request, fullpath=fullpath)
    messages.success(request, f"Copied file(s). Clipboard contains {len(copied)} files.")

    return redirect(reverse("data_view", kwargs=dict(uid=uid)))


@read_access(type=Job)
def job_file_copy(request, uid, path):
    # Get the root data where the file exists
    job = Job.objects.filter(uid=uid).first()
    fullpath = os.path.join(job.get_data_dir(), path)

    copied = auth.copy_file(request=request, fullpath=fullpath)
    messages.success(request, f"Copied file(s). Clipboard contains {len(copied)} files.")
    return redirect(reverse("job_view", kwargs=dict(uid=uid)))


@write_access(type=Project, fallback_view="recipe_list")
def recipe_paste(request, uid):
    """
    Pastes recipes from clipboard as a new recipes.
    """

    # The user performing the action.
    user = request.user

    # The project the paste will use.
    project = Project.objects.filter(uid=uid).first()

    # Contains the uids for the recipes that are to be copied.
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})

    paste_target = request.GET.get('target', const.COPIED_RECIPES)

    # Recipes in the clipboard are to be cloned
    paste_as_clone = paste_target == const.CLONED_RECIPES

    recipe_uids = clipboard.get(const.COPIED_RECIPES, [])

    # Select valid recipe uids.
    recipes = [Analysis.objects.filter(uid=uid).first() for uid in recipe_uids]

    # Keep existing recipes.
    recipes = filter(None, recipes)

    # The copy function for each recipe.
    def copy(instance):
        # Cascade the root if the recipe is being cloned
        if paste_as_clone:
            root = instance.root if instance.is_cloned else instance
        else:
            root = None
        recipe = auth.create_analysis(project=project, user=user, root=root,
                                      json_text=instance.json_text, security=instance.security,
                                      template=instance.template,
                                      name=instance.name, text=instance.text, stream=instance.image)
        return recipe

    # The list of new object created by the copy.
    new_recipes = list(map(copy, recipes))

    # Reset the session.
    clipboard[const.COPIED_RECIPES] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})

    # Notification after paste.
    messages.success(request, mark_safe(f"Pasted <b>{len(new_recipes)} recipes</b>  in clipboard"))

    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@write_access(type=Project, fallback_view="data_list")
def data_paste(request, uid):
    """Used to paste objects in results and data clipboards as a Data object."""
    project = Project.objects.filter(uid=uid).first()
    owner = request.user
    board = request.GET.get("board")
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    data_clipboard = clipboard.get(board, [])

    for datauid in data_clipboard:

        if board == const.COPIED_DATA:
            obj = Data.objects.filter(uid=datauid).first()
            dtype = obj.type
        else:
            obj = Job.objects.filter(uid=datauid).first()
            dtype = "DATA"

        if obj:
            auth.create_data(project=project, path=obj.get_data_dir(), user=owner, name=obj.name,
                             type=dtype, text=obj.text)

    clipboard[board] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})

    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@write_access(type=Project, fallback_view="data_list")
def file_paste(request, uid):
    project = Project.objects.filter(uid=uid).first()
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    file_clipboard = clipboard.get(const.COPIED_FILES, [])

    for single_file in file_clipboard:
        if os.path.exists(single_file):
            auth.create_data(project=project, path=single_file, user=request.user)

    clipboard[const.COPIED_FILES] = []
    request.session.update({settings.CLIPBOARD_NAME: clipboard})
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@read_access(type=Data)
def data_view(request, uid):
    "Show information specific to each data."

    data = Data.objects.filter(uid=uid).first()
    project = data.project

    context = dict(data=data, project=project, activate='Selected Data')
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

    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    form = forms.DataUploadForm(user=owner, project=project)

    if request.method == "POST":

        form = forms.DataUploadForm(data=request.POST, files=request.FILES, user=owner, project=project)

        if form.is_valid():
            data = form.save()
            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))

    uploaded_files = Data.objects.filter(owner=owner, method=Data.UPLOAD)

    # The current size of the existing data
    current_size = uploaded_files.aggregate(Sum("size"))["size__sum"] or 0

    # Maximum data that may be uploaded.
    maximum_size = owner.profile.max_upload_size * 1024 * 1024

    context = dict(project=project, form=form, activate="Add Data", maximum_size=maximum_size,
                   current_size=current_size)

    counts = get_counts(project)

    context.update(counts)

    return render(request, 'data_upload.html', context)



@read_access(type=Analysis)
def recipe_code_download(request, uid):
    """
    Download the raw recipe template as a file
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    try:
        # Fill in the script with json data.
        json_data = auth.fill_data_by_name(project=recipe.project, json_data=recipe.json_data)
        context = Context(json_data)
        script_template = Template(recipe.template)
        script = script_template.render(context)
    except Exception as exc:
        logger.error(exc)
        script = recipe.template

    # Trigger file download with name of the recipe
    filename = "_".join(recipe.name.split()) + ".sh"

    response = HttpResponse(script, content_type='text/plain')
    response['Content-Disposition'] = f'attachment; filename={filename}'

    return response


@read_access(type=Analysis)
@ratelimit(key='ip', rate='10/h', block=True, method=ratelimit.UNSAFE)
def recipe_run(request, uid):
    """
    View used to execute recipes and start a 'Queued' job.
    """

    analysis = Analysis.objects.filter(uid=uid).first()
    project = analysis.project

    # Form submission.
    if request.method == "POST":
        form = forms.RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data,
                                     data=request.POST, files=request.FILES)
        # The form validation will authorize the job.
        if form.is_valid():
            # Create the job from the recipe and incoming json data.
            job = auth.create_job(analysis=analysis, user=request.user, fill_with=form.cleaned_data)

            # Spool via UWSGI or start it synchronously.
            tasks.execute_job.spool(job_id=job.id)

            return redirect(reverse("job_view", kwargs=dict(uid=job.uid)))
    else:
        initial = dict(name=f"Results for: {analysis.name}")
        form = forms.RecipeInterface(request=request, analysis=analysis,
                                     json_data=analysis.json_data, initial=initial)

    is_runnable = auth.authorize_run(user=request.user, recipe=analysis)

    context = dict(project=project, analysis=analysis, form=form, is_runnable=is_runnable, activate='Run Recipe')
    context.update(get_counts(project))

    return render(request, 'recipe_run.html', context)


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

    return redirect(reverse('job_list', kwargs=dict(uid=job.project.uid)))


@read_access(type=Analysis)
def recipe_view(request, uid):
    """
    Edit meta-data associated with a recipe.
    """

    # The recipe that needs to be edited.
    recipe = Analysis.objects.filter(uid=uid).first()

    # The project that recipe belongs to.
    project = recipe.project
    user = request.user
    if request.method == "POST":
        # Form has been submitted
        form = forms.RecipeForm(data=request.POST, instance=recipe, files=request.FILES, user=user,
                                project=project)
        if form.is_valid():
            form.save()
            messages.success(request, "Editted Recipe")
            return redirect(reverse("recipe_view", kwargs=dict(uid=recipe.uid)))
    else:
        # Initial form loading via a GET request.
        form = forms.RecipeForm(instance=recipe, user=request.user, project=project)

    action_url = reverse('recipe_edit', kwargs=dict(uid=uid))
    initial = dict(name=f"Results for: {recipe.name}")

    run_form = forms.RecipeInterface(request=request, analysis=recipe,
                                     json_data=recipe.json_data, initial=initial)

    is_runnable = auth.authorize_run(user=request.user, recipe=recipe)

    recipe = Analysis.objects.filter(uid=uid).annotate(job_count=Count("job", filter=Q(job__deleted=False))).first()

    context = dict(recipe=recipe, project=project, form=form, is_runnable=is_runnable, name=recipe.name,
                   activate='Recipe View', run_form=run_form,
                   action_url=action_url)

    counts = get_counts(project)
    context.update(counts)
    return render(request, 'recipe_view.html', context)


@read_access(type=Project)
def recipe_create(request, uid):
    # Get the project

    project = Project.objects.filter(uid=uid).first()

    # Prepare the form
    name = "Recipe Name"
    initial = dict(name=name, uid=f'recipe-{util.get_uuid(5)}')
    form = forms.RecipeForm(user=request.user, initial=initial, project=project)

    if request.method == "POST":

        form = forms.RecipeForm(data=request.POST, creating=True, project=project, files=request.FILES,
                                user=request.user)
        if form.is_valid():

            image = form.cleaned_data['image']
            recipe_uid = form.cleaned_data['uid']
            name = form.cleaned_data['name']
            json_text = form.cleaned_data['json_text']
            template = form.cleaned_data['template']
            text = form.cleaned_data['text']
            rank = form.cleaned_data['rank']
            recipe = auth.create_analysis(uid=recipe_uid, stream=image, name=name, rank=rank,
                                          json_text=json_text, template=template,
                                          project=project, user=request.user, text=text)
            if request.user.is_superuser:
                security = Analysis.AUTHORIZED if form.cleaned_data.get('authorized') else Analysis.NOT_AUTHORIZED
                recipe.security = security
                recipe.save()

            return redirect(reverse("recipe_view", kwargs=dict(uid=recipe.uid)))

    action_url = reverse('recipe_create', kwargs=dict(uid=uid))
    context = dict(project=project, form=form, action_url=action_url,
                   activate='Create Recipe', name=name)
    counts = get_counts(project)
    context.update(counts)
    return render(request, 'recipe_view.html', context)


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
    clones = Analysis.objects.filter(root=recipe, deleted=False)

    if recipe.is_root and clones.exists():
        # Check if a root recipe
        msg = "Can not delete a cloned recipe."
        messages.success(request, msg)
    else:
        auth.delete_object(obj=recipe, request=request)
        msg = f"Deleted <b>{recipe.name}</b>." if recipe.deleted else f"Restored <b>{recipe.name}</b>."
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

    # The path is a GET parameter
    path = request.GET.get('path', "")

    stdout = job.stdout_log
    stderr = job.stderr_log
    if job.is_running():
        # Pass the most current stderr and stdout
        stdout_path = os.path.join(job.path, settings.JOB_STDOUT)
        stderr_path = os.path.join(job.path, settings.JOB_STDERR)
        stdout = open(stdout_path, 'r').read()
        stderr = open(stderr_path, 'r').read()

    context = dict(job=job, project=project, stderr=stderr, stdout=stdout,
                   activate='View Result', path=path)

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


def file_serve(request, path, obj):
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
    if size < 20:
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
    # Give user Share access.

    # Create access
    if access is None:
        Access.objects.create(user=user, project=project, access=Access.SHARE_ACCESS)

    # Update existing No Access to Share
    elif access.access == Access.NO_ACCESS:
        Access.objects.filter(id=access.id).update(access=Access.SHARE_ACCESS)

    messages.success(request, "Granted share access")
    return redirect(reverse('project_view', kwargs=dict(uid=project.uid)))


@login_required
def import_files(request, path=''):
    """
    Import files mounted on IMPORT_ROOT_DIR in settings
    """
    user = request.user

    if not user.profile.trusted:
        messages.error(request, 'Only trusted users may views this page.')
        return redirect(reverse('project_list'))

    current_path = os.path.abspath(os.path.join(settings.IMPORT_ROOT_DIR, path))

    if not current_path.startswith(settings.IMPORT_ROOT_DIR):
        messages.error(request, 'Outside root directory')
        rel_path = ''
    elif current_path == settings.IMPORT_ROOT_DIR:
        rel_path = ''
    else:
        rel_path = os.path.relpath(current_path, settings.IMPORT_ROOT_DIR)

    context = dict(rel_path=rel_path)

    return render(request, 'import_files.html', context=context)