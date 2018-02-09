from textwrap import dedent
import json
import logging
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.template import loader
from django.utils.safestring import mark_safe
from django.urls import reverse
from django.db.models import Q

from biostar import settings
from biostar.engine.models import Job, make_html, Project, Data, Analysis
from biostar.engine import auth


logger = logging.getLogger("engine")
register = template.Library()

JOB_COLORS = {
    Job.ZOMBIE: "orange", Job.SPOOLED: "violet",
    Job.ERROR: "red", Job.QUEUED: "teal",
    Job.RUNNING: "pink", Job.COMPLETED: "green"
}

DATA_COLORS = {
    Data.PENDING: "teal", Data.READY:"green", Data.ERROR :"red", Data.DELETED:""
}


@register.simple_tag
def sticky_label(obj):
    label = mark_safe('<span class ="ui label">Sticky</span>')
    return label if obj.sticky else ''


def access_denied_message(user, access):
    """
    Generates the access denied message
    """
    tmpl = loader.get_template('widgets/access_denied_message.html')
    context = dict(user=user, access=access)
    return tmpl.render(context=context)


@register.inclusion_tag('widgets/file_listing.html')
def file_listing(path, file_list, object, project_uid='', form=None):
    return dict(path=path, file_list=file_list, object=object, project_uid=project_uid, form=form)


def dir_url(object, path, current):

    url = object.url()
    path = path + "/" if path else ""

    if isinstance(object, Job):
        url = reverse("job_files_list", kwargs=dict(uid=object.uid, path=path + current.name))
    elif isinstance(object, Data):
        url = reverse("data_files_list", kwargs=dict(uid=object.uid, path=path + current.name))

    return mark_safe(f'<a href="{url}"><i class="folder icon"></i>{current.name}</a>')


@register.simple_tag
def input(path, current):
    path = path + "/" if path else ""
    return mark_safe(f'<input type="checkbox" name="paths" value="{path}{current.name}">')


@register.simple_tag
def file_url(object, path, current):

    if current.is_dir():
        return dir_url( object=object, path=path, current=current)

    media_url = settings.MEDIA_URL
    path = path + "/" if path else ""
    url = media_url + object.get_url() +  path + current.name

    return mark_safe(f'<a href="{url}"><i class="file text outline icon"></i>{current.name}</a>')


@register.filter
def has_data(request):
    """
    Checks if object in clipboard is a recipe.
    """

    uid = request.session.get("data_clipboard")
    data = Data.objects.filter(uid=uid).first()

    return bool(uid and data)


@register.filter
def has_files(request):
    "Checks if object in clipboard is a list of files belonging to a job"
    files = request.session.get("files_clipboard", [""])

    # Last item loaded into files clipboard is the job.uid the files belong to.
    job_uid = files.pop(-1)

    job = Job.objects.filter(uid=job_uid).first()
    if not job:
        return False

    # Some files in clipboard might be outside job path.
    files = [f for f in files if f.startswith(job.path)]

    # Files might still be empty at this point
    return True if files else False

def update_dict(iter):

    results=[dict(
                title=d.name,
                description=d.summary,
                url=d.url()) for d in iter
            ]

    return results

@register.inclusion_tag('widgets/search.html')
def search(request):
    #TODO: will probably be refractored after correctly using ajax
    #TODO: this is loaded in every page so there are a lot of queries going on even without a search,

    projects = auth.get_project_list(user=request.user)

    data = Data.objects.filter(~Q(state=Data.DELETED),project__in=projects)
    recipe = Analysis.objects.filter(~Q(state=Analysis.DELETED), project__in=projects)
    jobs = Job.objects.filter(~Q(state=Job.DELETED), project__in=projects)
    projects = projects.filter(~Q(state=Project.DELETED))

    content = [dict(name="Projects",results =update_dict(iter=projects)),
               dict(name="Data", results=update_dict(iter=data)),
               dict(name="Recipes", results=update_dict(iter=recipe)),
               dict(name="Jobs", results=update_dict(iter=jobs)),
               ]

    return dict(content=json.dumps(content))


@register.inclusion_tag('widgets/paste.html')
def paste(project, data=False, files=False):
    "Default provides template for pasting a recipe"

    action, url, message = "Paste Recipe", "recipe_paste", "a recipe"
    redir, board = "recipe_list", "recipe_clipboard"

    if data or files:
        # Change params for data or files

        action = "Paste Data" if data else "Paste File(s)"
        url = "data_paste" if data else "files_paste"
        board = "data_clipboard" if data else "files_clipboard"
        args, redir = dict(uid=project.uid), "data_list"
        message = "a data" if data else "files"

    paste_url = reverse(url, kwargs=dict(uid=project.uid))
    clear_url = reverse("clear_clipboard", kwargs=dict(uid=project.uid, redir=redir, board=board))

    return dict(action=action, paste_url=paste_url, clear_url=clear_url, message=message)


@register.filter
def has_recipe(request):
    """
    Checks if object in clipboard is a recipe.
    """
    uid = request.session.get("recipe_clipboard")
    recipe = Analysis.objects.filter(uid=uid).first()
    return bool(uid and recipe)


@register.inclusion_tag('widgets/pages.html')
def pages(instance):

    return dict(instance=instance)


@register.simple_tag
def privacy_label(project):
    label = mark_safe(f'<span class ="ui label">{project.get_privacy_display()}</span>' )
    return label


@register.inclusion_tag('widgets/authorization_required.html')
def security_label(analysis):
    context = dict(analysis=analysis)
    return context

@register.simple_tag
def job_color(job):
    """
    Returns a color based on job status.
    """
    return JOB_COLORS.get(job.state, "")

@register.simple_tag
def activate(value1, value2):
    """
    Returns a color based on job status.
    """
    return "active" if value1  == value2 else ''


@register.simple_tag
def data_color(data):
    "Return a color based on data status."

    return DATA_COLORS.get(data.state,  "")


@register.simple_tag
def type_label(data):
    if data.type:
        return mark_safe(f"<span class='ui label' > {data.type} </span>")
    return ""

@register.simple_tag
def state_label(data, error_only=False):

    label = f'<span class="ui { DATA_COLORS.get(data.state, "") } label"> {data.get_state_display()} </span>'

    # Error produce error only.
    if error_only and data.state != Data.ERROR:
        label = ""

    return mark_safe(label)


@register.simple_tag
def img(obj):
    """
    Returns the image associated with the object or a placeholder
    """
    if obj.image:
        return obj.image.url
    else:
        return static("images/placeholder.png")


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.inclusion_tag('widgets/project_name_bar.html')
def project_name_bar(project):
    """
    Returns a label for data sizes.
    """
    return dict(project=project)

@register.inclusion_tag('widgets/recipe_form.html')
def recipe_form(form):
    """
    Renders a recipe form.
    """
    return dict(form=form)

@register.inclusion_tag('widgets/access_form.html')
def access_form(project, user, form):
    """
    Generates an access form.
    """

    return dict(project=project, user=user, form=form)

@register.inclusion_tag('widgets/job_elapsed.html')
def job_minutes(job):
    """
    Returns a label for data sizes.
    """
    return dict(job=job)

@register.inclusion_tag('widgets/size_label.html')
def size_label(data):
    """
    Returns a label for data sizes.
    """
    return dict(data=data)


@register.inclusion_tag('widgets/form_errors.html')
def form_errors(form):
    """
    Turns form errors into a data structure
    """

    errorlist = [ ('', message) for message in form.non_field_errors() ]

    for field in form:
        for error in field.errors:
            errorlist.append((f'{field.name}:', error))

    context = dict(errorlist=errorlist)

    return context


@register.filter
def markdown(text):
    """
    Generates HTML from a markdown value.
    """
    if not text:
        return ''
    text = dedent(text)
    html = make_html(text)
    return mark_safe(html)


@register.inclusion_tag('widgets/menubar.html', takes_context=True)
def menubar(context, project=None, edit_project=False, create_project=False,
            data=None, edit_data=False, upload_data=False,
            analysis=None, edit_analysis=False, request=None
            ):
    user = context.request.user

    return dict(
        user=user,
        project=project, edit_project=edit_project, create_project=create_project,
        data=data, edit_data=edit_data, upload_data=upload_data,
        analysis=analysis, edit_analysis=edit_analysis, request=request
    )
