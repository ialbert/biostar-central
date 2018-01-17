from textwrap import dedent
import hjson, logging
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.template import loader
from django.forms import widgets
from django.db.models import Q
from django.utils.safestring import mark_safe
from django.urls import reverse

from biostar.engine import settings
from biostar.engine.models import Access, Job, make_html, Project, DataType, Data


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
def file_listing(path, file_list, object):

    return dict(path=path, file_list=file_list, object=object)


@register.simple_tag
def dir_url(object, path, current):

    url = object.url()
    path = path + "/" if path else ""

    if isinstance(object, Job):
        url = reverse("job_files_list", kwargs=dict(id=object.id, path=path + current.name))
    elif isinstance(object, Data):
        url = reverse("data_files_list", kwargs=dict(id=object.id, path=path + current.name))

    return mark_safe(f'<a href="{url}"><i class="folder icon"></i>{current.name}</a>')


@register.simple_tag
def file_url(object, path, current):

    media_url = settings.MEDIA_URL
    path = path + "/" if path else ""
    url = media_url + object.get_url() +  path + current.name

    return mark_safe(f'<a href="{url}"><i class="file text outline icon"></i>{current.name}</a>')



@register.inclusion_tag('widgets/copy_interface.html')
def copy_interface(form, projects, duplicate=False):
    """Copy an instance from (from_id) to a list of allowed projects"""

    copying = "recipe" if duplicate else "data"
    return dict(projects=projects, form=form, duplicate=duplicate, copying=copying)

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
def access_label(project, user):

    if user.is_anonymous:
        return ""

    access = Access.objects.filter(project=project, user=user).first()

    if not access and project.privacy == Project.PUBLIC:
        return ""

    label = mark_safe(f'<span class ="ui green label">{access.get_access_display()}</span>')
    return label


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


@register.filter
def single_file(data):
    "Return true if data only contains one file"

    if len(data.get_files()) > 1:
        return False

    return True


@register.inclusion_tag('widgets/project_name_bar.html')
def project_name_bar(project):
    """
    Returns a label for data sizes.
    """
    return dict(project=project)

@register.inclusion_tag('interface/recipe_form.html')
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


@register.inclusion_tag('widgets/type_label.html')
def type_label(data):
    """
    Returns a label for a data type.
    """

    color = ""
    query = DataType.objects.filter(project=data.project, symbol=data.data_type).first()

    if not query:
        query, color = f"{data.data_type} Not Recognized", "yellow"

    return dict(label=query, color=color)


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

    #if errorlist:
    #    print (errorlist, type(errorlist[0]), dir(errorlist[0]))

    return context


@register.simple_tag
def field_state(field):
    """
    Returns the error label for a field.
    """

    if field.errors:
        return 'error'
    else:
        return ''


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


@register.inclusion_tag('widgets/breadcrumb.html')
def breadcrumb(steps):
    """
    Generates the breadcrumb for a page.
    """
    return dict(steps=steps)


@register.inclusion_tag('widgets/menubar.html', takes_context=True)
def menubar(context, project=None, edit_project=False, create_project=False,
            data=None, edit_data=False, upload_data=False,
            analysis=None, edit_analysis=False
            ):
    user = context.request.user

    return dict(
        user=user,
        project=project, edit_project=edit_project, create_project=create_project,
        data=data, edit_data=edit_data, upload_data=upload_data,
        analysis=analysis, edit_analysis=edit_analysis,
    )
