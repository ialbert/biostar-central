import glob
import logging
import json

import mistune
from django.conf import settings
from django.contrib.auth.decorators import user_passes_test
from django.core.paginator import Paginator
from django.shortcuts import render, redirect
from django.urls import reverse
from sendfile import sendfile

# from django.utils.safestring import mark_safe
from biostar.breadcrumb import breadcrumb_builder
from . import tasks
from .decorators import object_access
from .forms import *
from .models import (Project, Data, Analysis, Job, Access)


def join(*args):
    return os.path.abspath(os.path.join(*args))


# Objects per page when looking at lists
OBJ_PER_PAGE = 10

# The current directory
__CURRENT_DIR = os.path.dirname(__file__)
__DOCS_DIR = join(__CURRENT_DIR, "docs")


def valid_path(path):
    path = os.path.abspath(path)
    return path.startswith(__DOCS_DIR)


logger = logging.getLogger('engine')


def make_html(text):
    return mistune.markdown(text)


def pages(request, instance):
    paginator = Paginator(instance, OBJ_PER_PAGE)
    page = request.GET.get('page', 1)

    return paginator.page(page)


def docs(request, name):
    patt = join(__DOCS_DIR, name) + ".*"
    files = glob.glob(patt)
    if not files:
        msg = f"Cannot be find the requested page: {name} "
        messages.error(request, msg)
        return redirect("index")
    if len(files) > 1:
        msg = f"Multiple files match: {{name}}"
        messages.warning(request, msg)
    target = files[0]
    content = open(target).read()

    # Render markdown into HTML.
    if target.endswith(".md"):
        content = make_html(content)

    title = name.replace("-", " ").replace("_", " ").title()
    context = dict(content=content, title=title, steps=[])
    return render(request, 'info/doc_base.html', context=context)


def index(request):
    steps = breadcrumb_builder([HOME_ICON])
    context = dict(steps=steps)
    return render(request, 'index.html', context)


@user_passes_test(lambda u: u.is_superuser)
def site_admin(request):
    '''
    Administrative view. Lists the admin project and job.
    '''
    steps = breadcrumb_builder([HOME_ICON])
    projects = Project.objects.all()
    context = dict(steps=steps, projects=projects)
    return render(request, 'admin_index.html', context=context)



def clear_clipboard(request, uid, redir="project_view"):
    "Clear copy object held in clipboard"

    request.session["clipboard"] = None

    return redirect(reverse(redir, kwargs=dict(uid=uid)))

@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def project_users(request, uid):
    """
    Manage project users
    """

    # project = Project.objects.filter(uid=uid).first()
    #
    # # Search query
    # q = request.GET.get("q")
    #
    # # Users already with access to current project
    # users = [access.user for access in project.access_set.all() if access.access > Access.NO_ACCESS]
    #
    # # Users that have been searched for.
    # targets = []
    #
    # steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ADD_USER],
    #                            project=project)
    #
    # if request.method == "POST":
    #     form = ChangeUserAccess(data=request.POST)
    #     if form.is_valid():
    #         form.change_access()
    #         messages.success(request, "Changed access to this project")
    #         return redirect(reverse("project_users", kwargs=dict(uid=project.uid)))
    #
    #     messages.error(request, mark_safe(form.non_field_errors()))
    #
    # if q:
    #     targets = User.objects.filter(Q(email__contains=q) | Q(first_name__contains=q))
    #
    # current = access_forms(users=users, project=project)
    # results = access_forms(users=targets, project=project)
    # context = dict(steps=steps, current=current, project=project, results=results)
    # return render(request, "project_users.html", context=context)
    return redirect(reverse("project_view", kwargs=dict(uid=uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def project_types(request, uid):
    "Manage data types belonging to a project from a project"

    # project = Project.objects.filter(uid=uid).first()
    # steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, PROJECT_TYPES],
    #                            project=project)
    # form = CreateDataTypeForm(project=project)
    #
    # if request.method == "POST":
    #     form = CreateDataTypeForm(project=project, data=request.POST)
    #     if form.is_valid():
    #         form.save()
    #         return redirect(reverse("project_types", kwargs=dict(uid=project.uid)))
    #
    #     messages.error(request, mark_safe(form.errors))
    #
    # current = project.datatype_set.all()
    # context = dict(project=project, form=form, steps=steps, current=current)
    # return render(request, "project_types.html", context=context)
    return redirect(reverse("project_view", kwargs=dict(uid=uid)))


def project_list(request):
    projects = auth.get_project_list(user=request.user).order_by("-sticky", "-privacy")
    projects = projects.order_by("-privacy", "-sticky", "-date", "-id")
    projects = pages(request, instance=projects)

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])
    context = dict(projects=projects, steps=steps)

    return render(request, "project_list.html", context)


def data_list(request, uid):
    """
    Returns the list of data for a project id.
    """

    return project_view(request=request, uid=uid, template_name="data_list.html", active='data')


def recipe_list(request, uid):
    """
    Returns the list of recipes for a project id.
    """
    return project_view(request=request, uid=uid, template_name="recipe_list.html", active='recipes')


def job_list(request, uid):
    """
    Returns the list of recipes for a project id.
    """
    return project_view(request=request, uid=uid, template_name="job_list.html", active='jobs')


def files_list(request, instance, steps, error_redirect, template_name, path='',
               extra_context={}):
    "File navigator used for jobs and data"

    # Instance is expected to be a Job or Data object.
    exclude = ''
    if isinstance(instance, Job):
        root = instance.path
    else:
        # Exclude toc from file_list
        exclude = os.path.basename(instance.get_path())
        root = instance.get_data_dir()

    target_path = join(root, path)

    if not target_path.startswith(root) or (not os.path.exists(target_path)):
        # Attempting to access a file outside of the job directory
        messages.error(request, "Path not in directory.")
        return redirect(error_redirect)

    #TODO: can not exclude toc from other projects when copying
    # These are pathlike objects with attributes such as name, is_file
    file_list = list(filter(lambda p: p.name != exclude, os.scandir(target_path)))

    # Sort by properties
    file_list = sorted(file_list, key=lambda p: (p.is_file(), p.name))
    context = dict(file_list=file_list, instance=instance, steps=steps, path=path)
    context.update(extra_context)

    return render(request, template_name, context)


def get_counts(project):
    data_count = Data.objects.filter(project=project).count()
    recipe_count = Analysis.objects.filter(project=project).count()
    result_count = Job.objects.filter(project=project).count()
    return dict(
        data_count=data_count, recipe_count=recipe_count, result_count=result_count
    )


@object_access(type=Project, access=Access.READ_ACCESS)
def project_view(request, uid, template_name="recipe_list.html", active='recipes'):
    user = request.user

    project = Project.objects.filter(uid=uid).first()

    # Project not found.
    if not project:
        messages.error(request, "Project not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON],
                               project=project)

    # Show counts for the project.
    counts = get_counts(project)

    # Select all the data in the project.
    data_list = Data.objects.filter(project=project).order_by("sticky", "-date").all()
    recipe_list = Analysis.objects.filter(project=project).order_by("-date").all()
    job_list = Job.objects.filter(project=project).order_by("-date").all()

    # Filter job results by analysis
    filter = request.GET.get('filter', '')
    if filter:
        filter = Analysis.objects.filter(id=filter).first()
        job_list = job_list.filter(analysis=filter)

    if user.is_authenticated():
        access = Access.objects.filter(user=user, project=project).first()
    else:
        access = None

    context = dict(project=project, access=access, data_list=data_list, recipe_list=recipe_list, job_list=job_list,
                   active=active, steps=steps, filter=filter, help_text=project.html)
    context.update(counts)

    return render(request, template_name, context)


@object_access(type=Project, access=Access.EDIT_ACCESS, url='project_view')
def project_edit(request, uid):
    project = auth.get_project_list(user=request.user).filter(uid=uid).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON], project=project)
    form = ProjectForm(instance=project)

    if request.method == "POST":
        form = ProjectForm(request.POST, request.FILES, instance=project)
        if form.is_valid():
            form.save()
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

        messages.error(request, mark_safe(form.errors))

    context = dict(project=project, steps=steps, form=form)
    return render(request, "project_edit.html", context=context)


def project_create(request):
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON])

    if request.user.is_anonymous:
        messages.error(request, "You must be logged in to create a project.")
        return redirect(reverse("project_list"))

    initial = dict(name="Project Name", text="project description", summary="project summary")
    form = ProjectForm(initial=initial)

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
            uid = form.cleaned_data["uid"]
            owner = request.user
            project = auth.create_project(user=owner, name=name, summary=summary, text=text,
                                          stream=stream, sticky=sticky, privacy=privacy,
                                          uid=uid)
            project.save()
            return redirect(reverse("project_view", kwargs=dict(uid=project.uid)))

        messages.error(request, mark_safe(form.errors))

    context = dict(steps=steps, form=form)
    return render(request, "project_create.html", context=context)

@object_access(type=Data, access=Access.READ_ACCESS)
def data_view(request, id):
    data = Data.objects.filter(id=id).first()

    if not data:
        messages.error(request, "Data not found.")
        logger.error(f"data.id={id} looked for but not found.")
        return redirect(reverse("project_list"))

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_ICON],
                               project=data.project, data=data)

    project = data.project
    context = dict(data=data, steps=steps, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "data_view.html", context)


@object_access(type=Data, access=Access.READ_ACCESS, url='data_view')
def data_copy(request, id):
    "Store Data object in request.sessions['clipboard'] "


    data = Data.objects.filter(pk=id).first()
    project = data.project

    request.session["clipboard"] = data.uid
    messages.success(request, f"Copied {data.name} to Clipboard")

    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def data_paste(request, uid):
    "Paste data stored in the clipboard into project"

    data_uid = request.session.get("clipboard")
    data = Data.objects.filter(uid=data_uid).first()
    project = Project.objects.filter(uid=uid).first()

    if not data:
        messages.error(request, "Data not found.")
        return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))

    data_files = data.get_files()

    if len(data_files) > 1:
        # Link the dir if more than one file present
        data_files = [ join(data_files[0], "..") ]

    path = data_files.pop()
    # Create new data by linking file(s)
    new_data = auth.create_data(project=project, name=f"Copy of {data.name}", path=path,
                                summary=data.summary, text=data.text,
                                data_type=data.data_type)
    new_data.save()

    # Clear clipboard
    request.session["clipboard"] = None

    messages.success(request, f"Pasted {data.name} to {project.name}")
    return redirect(reverse("data_list", kwargs=dict(uid=project.uid)))


@object_access(type=Data, access=Access.EDIT_ACCESS, url='data_view')
def data_edit(request, id):
    data = Data.objects.filter(id=id).first()
    project = data.project
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)
    form = DataEditForm(instance=data, project=project, initial=dict(data_type=data.data_type))

    if request.method == "POST":
        form = DataEditForm(data=request.POST, instance=data, project=project)
        if form.is_valid():
            form.save()
            return redirect(reverse("data_view", kwargs=dict(id=data.id)))

        messages.error(request, mark_safe(form.errors))

    context = dict(data=data, steps=steps, form=form)
    return render(request, 'data_edit.html', context)


@object_access(type=Project, access=Access.UPLOAD_ACCESS, url='data_list')
def data_upload(request, uid):
    owner = request.user
    project = Project.objects.filter(uid=uid).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_UPLOAD],
                               project=project)
    form = DataUploadForm(project=project)
    if request.method == "POST":
        form = DataUploadForm(data=request.POST, files=request.FILES, project=project)

        if form.is_valid():
            text = form.cleaned_data["text"]
            stream = form.cleaned_data["file"]
            summary = form.cleaned_data["summary"]
            data_type = form.cleaned_data["data_type"]
            name = stream.name
            data = auth.create_data(stream=stream, name=name,
                                    text=text, user=owner, project=project, summary=summary,
                                    data_type=data_type)
            messages.info(request, f"Uploaded: {data.name}. Edit the data to set its type.")
            return redirect(reverse("data_list", kwargs={'uid': project.uid}))

        messages.error(request, mark_safe(form.errors))

    context = dict(project=project, steps=steps, form=form)
    return render(request, 'data_upload.html', context)


@object_access(type=Data, access=Access.ADMIN_ACCESS, url='data_view')
def data_files_list(request, id, path=''):

    data = Data.objects.filter(id=id).first()
    project = data.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, DATA_LIST_ICON, DATA_ICON],
                               project=project, data=data)

    context = dict(activate='selection', data=data, project=project)

    counts = get_counts(project)
    context.update(counts)

    return files_list(request=request, instance=data, path=path,steps=steps,
                      error_redirect=reverse("data_files_entry", kwargs=dict(id=id)),
                      template_name="data_files_list.html", extra_context=context)


@object_access(type=Data, access=Access.ADMIN_ACCESS, url='data_view')
def data_download(request, id):
    "Download data found in a project"

    data = Data.objects.filter(id=id).first()

    if not data:
        messages.error(request, "Data Not Found")
        return redirect(reverse("project_list"))

    data_file = data.get_files()

    if len(data_file) > 1:
        # Redirect to file navigator if there are multiple files
        return redirect(reverse("data_files_entry", kwargs=dict(id=data.id)))

    # Only one file expected at this point
    data_file = data_file.pop()

    if not os.path.exists(data_file):
        messages.error(request, "Data object does not contain a valid file")
        return redirect(reverse("data_view", kwargs=dict(id=id)))

    return sendfile(request, data_file, attachment=f"{data.name}")


@object_access(type=Analysis, access=Access.READ_ACCESS)
def recipe_view(request, id):
    """
    Returns an analysis view based on its id.
    """
    analysis = Analysis.objects.filter(id=id).first()
    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, ANALYSIS_LIST_ICON,
                                ANALYSIS_VIEW_ICON], project=analysis.project, analysis=analysis)

    project = analysis.project

    context = dict(analysis=analysis, steps=steps, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "recipe_view.html", context)


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
def recipe_run(request, id):
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                ANALYSIS_VIEW_ICON, ANALYSIS_RUN_ICON],
                               project=project, analysis=analysis)

    if request.method == "POST":
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, data=request.POST)

        if form.is_valid():

            # The desired name of for the results.
            name = form.cleaned_data.get("name")

            # Generates the JSON data from the bound form field.
            json_data = form.fill_json_data()

            # Create the job from the json.
            state = Job.SPOOLED if tasks.HAS_UWSGI else Job.QUEUED
            job = auth.create_job(analysis=analysis, user=request.user, json_data=json_data, name=name, state=state)

            # Spool the job right away if UWSGI exists.
            if tasks.HAS_UWSGI:
                tasks.execute_job.spool(job_id=job.id)

            return redirect(reverse("job_list", kwargs=dict(uid=project.uid)))
    else:
        initial = dict(name=f"Results for: {analysis.name}")
        form = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=initial)

    context = dict(project=project, analysis=analysis, steps=steps, form=form, activate='selection')
    counts = get_counts(project)
    context.update(counts)

    return render(request, 'recipe_run.html', context)


@object_access(type=Analysis, access=Access.READ_ACCESS, url='recipe_view')
def recipe_copy(request, id):
    "Store Analysis object in request.sessions['clipboard'] "

    recipe = Analysis.objects.filter(pk=id).first()
    project = recipe.project

    request.session["clipboard"] = recipe.uid
    messages.success(request, f"Copied {recipe.name} to clipboard")

    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@object_access(type=Project, access=Access.ADMIN_ACCESS, url='project_view')
def recipe_paste(request, uid):
    "Paste recipe stored in the clipboard into project"

    recipe_uid = request.session.get("clipboard")
    recipe = Analysis.objects.filter(uid=recipe_uid).first()
    project = Project.objects.filter(uid=uid).first()

    if not recipe:
        messages.error(request, "Recipe not found.")
        return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))

    attrs = auth.get_analysis_attr(recipe, project=project)
    attrs.update(stream=recipe.image, name=f"Copy of {recipe.name}", security=recipe.security)
    new_recipe = auth.create_analysis(**attrs)
    new_recipe.save()

    messages.success(request, f"Pasted {recipe.name} to {project.name}")
    request.session["clipboard"] = None

    return redirect(reverse("recipe_list", kwargs=dict(uid=project.uid)))


@object_access(type=Analysis, access=Access.RECIPE_ACCESS, url='recipe_view')
def recipe_code(request, id):
    """
    Displays and allows edit on a recipe code.

    Because we allow a preview even for unauthenicated users the view
    is a lot more complicated than a typical DJANO form handler.
    """
    user = request.user

    # There has to be a recipe to work with.
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    name = analysis.name

    # This is the navbat.
    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

    if request.method == "POST":
        form = EditCode(user=user, project=project, data=request.POST)

        if form.is_valid():

            fixend = lambda t: t.replace('\r\n', '\n')

            template = fixend(form.cleaned_data['template'])

            # Preview action will let the form cascade through.
            save = form.cleaned_data['action'] == 'SAVE'

            # The changes will commited on SAVE only.

            analysis.json_text = fixend(form.cleaned_data['json'])

            # Changes to template will require a review ( only when saving ).
            if auth.template_changed(analysis=analysis, template=template) and save:
                # Switch on the untrusted flag when the template changes.
                analysis.security = Analysis.UNDER_REVIEW

            # Set the new template.
            analysis.template = template

            # The SAVE action commits the changes on the analysis.
            if save:
                analysis.save()
                messages.info(request, "The recipe code has been updated.")
                return redirect(reverse("recipe_view", kwargs=dict(id=analysis.id)))
    else:
        # This gets triggered on a GET request.
        initial = dict(template=analysis.template, json=analysis.json_text)
        form = EditCode(user=user, project=project, initial=initial)

    # Bind the JSON to the form.
    recipe = RecipeInterface(request=request, analysis=analysis, json_data=analysis.json_data, initial=dict(name=name))

    # This generates a "fake" unsaved job.
    job = auth.create_job(analysis=analysis, json_data=analysis.json_data, save=False)

    # Create the script for the "fake" job.
    data, script = auth.generate_script(job)

    # Populate the context.
    context = dict(project=project, analysis=analysis, steps=steps, form=form, script=script, recipe=recipe)
    return render(request, 'recipe_code.html', context)


@object_access(type=Project, access=Access.EDIT_ACCESS, url='recipe_list')
def recipe_create(request, uid):
    """
    Create recipe with empty template and json spec
    """

    # project = Project.objects.filter(uid=uid).first()
    #
    # steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON], project=project)
    # action_url = reverse('recipe_create', kwargs=dict(uid=project.uid))
    # back_url = reverse('recipe_list', kwargs=dict(uid=project.uid))
    #
    # if request.method == "POST":
    #     form = RecipeForm(data=request.POST, files=request.FILES)
    #
    #     if form.is_valid():
    #         # Empty Analysis Template is authorized on creation
    #         security = Analysis.AUTHORIZED
    #         name = form.cleaned_data["name"]
    #         text = form.cleaned_data["text"]
    #         summary = form.cleaned_data["summary"]
    #         stream = form.cleaned_data["image"]
    #         sticky = form.cleaned_data["sticky"]
    #
    #         recipe = auth.create_analysis(project=project, json_text="{}", template="",
    #                                       user=request.user, summary=summary, name=name, text=text,
    #                                       security=security, stream=stream, sticky=sticky)
    #         recipe.save()
    #         messages.success(request, "Recipe created")
    #         return redirect(back_url)
    #
    # form = RecipeForm()
    # context = dict(steps=steps, analysis={"name": "New Analysis"},
    #                project=project, form=form, action_url=action_url, back_url=back_url)
    # return render(request, 'recipe_edit.html', context)
    return redirect(reverse("project_list", kwargs=dict(uid=uid)))


@object_access(type=Analysis, access=Access.EDIT_ACCESS, url='recipe_view')
def recipe_edit(request, id):
    "Edit recipe Info"
    analysis = Analysis.objects.filter(id=id).first()
    project = analysis.project

    steps = breadcrumb_builder([PROJECT_ICON, ANALYSIS_LIST_ICON, ANALYSIS_VIEW_ICON,
                                ANALYSIS_RECIPE_ICON], project=project, analysis=analysis)

    action_url = reverse('recipe_edit', kwargs=dict(id=analysis.id))
    back_url = reverse('recipe_view', kwargs=dict(id=analysis.id))

    if request.method == "POST":
        form = RecipeForm(data=request.POST, files=request.FILES, instance=analysis)
        if form.is_valid():
            recipe = form.save()
            return redirect(reverse("recipe_view", kwargs=dict(id=recipe.id)))

        messages.error(request, mark_safe(form.errors))

    form = RecipeForm(instance=analysis)
    context = dict(steps=steps, analysis=analysis, project=project, form=form, action_url=action_url, back_url=back_url)

    return render(request, 'recipe_edit.html', context)


@object_access(type=Job, access=Access.EDIT_ACCESS, url="job_view")
def job_edit(request, id):
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON, RESULT_LIST_ICON,
                                RESULT_VIEW_ICON], job=job, project=project)

    if request.method == "POST":
        form = JobEditForm(data=request.POST, files=request.FILES, instance=job)
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

    steps = breadcrumb_builder([HOME_ICON, PROJECT_LIST_ICON, PROJECT_ICON,
                                RESULT_VIEW_ICON], job=job, project=project)

    project = job.project

    context = dict(job=job, steps=steps, project=project, activate='selection')

    counts = get_counts(project)
    context.update(counts)

    return render(request, "job_view.html", context=context)


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_result_view(request, id):
    """
    Returns the primary result of a job.
    """

    job = Job.objects.filter(id=id).first()
    index = job.json_data.get("settings", {}).get("index")

    if job.state == Job.COMPLETED:
        url = reverse("job_files_entry", kwargs=dict(id=id))

        if index:
            url = settings.MEDIA_URL + job.get_url(path=index)

        return redirect(url)

    return redirect(reverse("job_view", kwargs=dict(id=id)))


def block_media_url(request, **kwargs):
    "Block users from accessing media directory using urls"

    messages.error(request, f"Not allowed")
    return redirect(reverse("project_list"))


@object_access(type=Job, access=Access.READ_ACCESS, url="job_view")
def job_files_list(request, id, path=''):
    """
    Returns the directory view of the job.
    """
    job = Job.objects.filter(id=id).first()
    project = job.project

    steps = breadcrumb_builder(
        [PROJECT_LIST_ICON, PROJECT_ICON, RESULT_VIEW_ICON, RESULT_INDEX_ICON],
        job=job, project=project)

    context = dict(activate='selection', job=job, project=project)

    counts = get_counts(project)
    context.update(counts)

    return files_list(request=request, instance=job, path=path, steps=steps,
                      template_name="job_files_list.html",
                      error_redirect=reverse("job_files_entry", kwargs=dict(id=id)),
                      extra_context=context)

