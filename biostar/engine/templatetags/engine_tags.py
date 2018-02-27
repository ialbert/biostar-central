from textwrap import dedent
import json
import logging
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.template import defaultfilters
from django.utils.safestring import mark_safe
from django.urls import reverse
from django.db.models import Q
from django.test.client import RequestFactory

from biostar import settings
from biostar.engine.models import Job, make_html, Project, Data, Analysis, Access
from biostar.engine import auth


logger = logging.getLogger("engine")
register = template.Library()

JOB_COLORS = {
    Job.SPOOLED: "violet",
    Job.ERROR: "red", Job.QUEUED: "teal",
    Job.RUNNING: "orange", Job.COMPLETED: "green"
}

DATA_COLORS = {
    Data.PENDING: "teal", Data.READY:"green", Data.ERROR :"red"
}


@register.simple_tag
def sticky_label(obj):
    label = mark_safe('<span class ="ui label">Sticky</span>')
    return label if obj.sticky else ''



@register.inclusion_tag('widgets/file_listing.html')
def file_listing(path, file_list, object, project_uid='', form=None):
    return dict(path=path, file_list=file_list, object=object, project_uid=project_uid, form=form)


@register.inclusion_tag('widgets/action_bar.html', takes_context=True)
def action_bar(context, instance, edit_url, extra=""):

    edit_url = reverse(edit_url, kwargs=dict(uid=instance.uid))

    if isinstance(instance, Job):
        obj_type = "job"
    elif isinstance(instance, Data):
        obj_type = "data"
    elif isinstance(instance, Analysis):
        obj_type = "recipe"
    else:
        obj_type = None

    action_url = reverse("toggle_state", kwargs=dict(uid=instance.uid, obj_type=obj_type))
    return dict(instance=instance, edit_url=edit_url, user=context["user"],
                action_url=action_url, extra=extra)


def dir_url(object, path, current):

    url = object.url()
    path = path + "/" if path else ""

    if isinstance(object, Job):
        url = reverse("job_browser", kwargs=dict(uid=object.uid, path=path + current.name))
    elif isinstance(object, Data):
        url = reverse("data_browser", kwargs=dict(uid=object.uid, path=path + current.name))

    return mark_safe(f'<a href="{url}"><i class="folder icon"></i>{current.name}</a>')


@register.inclusion_tag('widgets/list_view.html')
def list_view(projects=None, data_list=None, recipe_list=None, job_list=None):

    return dict(projects=projects, data_list=data_list, recipe_list=recipe_list,
                job_list=job_list)

@register.simple_tag
def input(path, current):

    path = path + "/" if path else ""
    input_template = f'<input type="checkbox" name="paths" value="{path}{current.name}">'
    return mark_safe(input_template)


@register.simple_tag
def file_url(object, path, current):
    """
    Used in file navigator to render file info in template.
    param: current is
    """

    assert isinstance(object, Data) or isinstance(object, Job)
    if current.is_dir():
        return dir_url( object=object, path=path, current=current)

    path = path + "/" if path else ""
    if isinstance(object, Data):
        url = reverse("data_serve", kwargs=dict(uid=object.uid, path=path + current.name))
    else:
        url = reverse("job_serve", kwargs=dict(uid=object.uid, path=path + current.name))

    byte = current.stat(follow_symlinks=True).st_size

    # Make the size user friendly
    size = f"{defaultfilters.filesizeformat(byte)}"
    file_template = f"""
            <a href="{url}"><i class="file text outline icon"></i>{current.name}
            <span class="ui right floated mini label">{size}</span>
            </a>
            """
    return mark_safe(file_template)


@register.filter
def has_data(request):
    """
    Checks if object in clipboard is a recipe.
    """

    uid = request.session.get("data_clipboard")

    # Dump data clipboard without resting it
    data_list = auth.dump_data_clipboard(request=request)

    return bool(uid and data_list)

@register.filter
def is_checkbox(field):
    "Check if current field is a checkbox"

    return True if field.field.widget.input_type == "checkbox" else False


@register.filter
def has_files(request):
    "Checks if object in clipboard is a list of files belonging to a job"
    files = request.session.get("files_clipboard", None)

    # Last item loaded into files clipboard is the instance.uid the files belong to.
    root_path = auth.validate_files_clipboard(request=request, clipboard=files)
    if not root_path:
        return False

    # Some files in clipboard might be outside job path.
    files = [f for f in files if f.startswith(root_path)]

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

    projects = auth.get_project_list(user=request.user)

    data = Data.objects.filter(deleted=False,project__in=projects)
    recipe = Analysis.objects.filter(deleted=False, project__in=projects)
    jobs = Job.objects.filter(deleted=False, project__in=projects)
    projects = projects.filter(deleted=False)

    content = [dict(name="Projects",results =update_dict(iter=projects)),
               dict(name="Data", results=update_dict(iter=data)),
               dict(name="Recipes", results=update_dict(iter=recipe)),
               dict(name="Jobs", results=update_dict(iter=jobs)),
               ]

    return dict(content=json.dumps(content))


@register.inclusion_tag('widgets/paste.html')
def paste(project, request=None, board=''):
    "Default provides template for pasting a recipe"

    files = board == 'files_clipboard'
    active_board = [] if not request else request.session.get(board, [])
    # Last time in the files_clipboard is an instance.uid, not a file.
    ndata = len(active_board) - 1 if files and active_board else len(active_board)

    # Map a clipboard to respective information
    info = dict(
        recipe_clipboard={'url':"recipe_paste", 'nactive':"a recipe", 'redirl':"recipe_list"},
        data_clipboard={'url':"data_paste", 'nactive':f"{ndata} data",'redir':"data_list"},
        files_clipboard={'url':"files_paste", 'nactive':f"{ndata} file(s)", 'redir':"data_list"}
    )

    view_name = info.get(board, {}).get('url')
    paste_url = reverse(view_name, kwargs=dict(uid=project.uid))
    clear_url = reverse("clear_clipboard", kwargs=dict(uid=project.uid,
                                                       url=info.get(board, {}).get('redir'),
                                                       board=board))

    return dict(paste_url=paste_url, clear_url=clear_url, message=info.get(board, {}).get('nactive'))


@register.filter
def has_recipe(request):
    """
    Checks if object in clipboard is a recipe.
    """
    uid = request.session.get("recipe_clipboard")
    recipe = Analysis.objects.filter(uid=uid).first()
    return bool(uid and recipe)



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
def access_color(user, project):

    if user.is_authenticated:
        access = Access.objects.filter(user=user, project=project).first()
    else:
        access = None

    if access and access.access in (Access.WRITE_ACCESS, Access.OWNER_ACCESS):
        return "green"
    else:
        return ""


@register.simple_tag
def type_label(data):

    if data.type:
        label = lambda x:f"<span class='ui label' > {x} </span>"
        types = [label(t) for t in data.type.split(',')]
        return mark_safe(''.join(types))
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


@register.inclusion_tag('widgets/project_name_bar.html', takes_context=True)
def project_name_bar(context, project):
    """
    Returns a label for project.
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


@register.inclusion_tag('widgets/project_action_bar.html')
def project_action_bar(user, project):
    return dict(use=user, project=project)

@register.inclusion_tag('widgets/job_elapsed.html')
def job_minutes(job):

    return dict(job=job)


@register.simple_tag
def size_label(data):
    """
    Returns a label for data sizes.
    """
    size = f"{defaultfilters.filesizeformat(data.size)}"
    return mark_safe(f"<span class='ui mini label'>{size}</span>")


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
def menubar(context, request=None):
    user = context.request.user

    return dict(user=user, request=request)
