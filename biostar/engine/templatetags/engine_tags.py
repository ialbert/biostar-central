from textwrap import dedent
import hjson, logging
from django import forms
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.template import Template, Context
from django.forms import widgets

from django.utils.safestring import mark_safe

from biostar.engine import const
from biostar.engine.models import Project, Job, make_html
from biostar.engine import factory

logger = logging.getLogger("engine")
register = template.Library()

JOB_COLORS = {
    Job.ZOMBIE: "orange", Job.SPOOLED: "pink",
    Job.ERROR: "red", Job.QUEUED: "blue", Job.RUNNING: "teal", Job.COMPLETED: "green"
}


def make_form_field(data, project=None):

    display_type = data.get("display_type")

    # Fields with no display type are not visible.
    if not display_type:
        return

    # Uploaded data is accessed via paths or links.
    path_or_link = data.get("path") or data.get("link")

    if path_or_link and project:
        # Project specific data needs a special field.
        data_type = data.get("data_type")
        field = factory.data_field_generator(data, project=project, data_type=data_type)
    else:

        func = factory.TYPE2FUNC.get(display_type)
        if not func:
            logger.error(f"Invalid display_type={display_type}")
            return
        field = func(data)

    return field


@register.inclusion_tag('widgets/json_form.html')
def generate_fields(json_text, project=None, form=None):

    fields = []

    json_data = hjson.loads(json_text)

    for name, data in json_data.items():
        field = make_form_field(data, project)
        if field:
            field.widget.attrs["name"] = name
            # Returns <django.forms.fields.CharField object> instead of html if the field isnt
            # bound to a form
            if form:
                field = forms.forms.BoundField(form=form, field=field, name=name)
            else:
                field = {"field": field}
            fields.append(field)

    return dict(fields=fields)


@register.filter
def generate_script(template, json_text):

    json_data = hjson.loads(json_text)
    template = Template(template)
    context = Context(json_data)

    return template.render(context)


@register.simple_tag
def sticky_label(obj):
    label = mark_safe('<span class ="ui label">Sticky</span>')
    return label if obj.sticky else ''



@register.simple_tag
def privacy_label(project):
    label = mark_safe(f'<span class ="ui label">{project.get_privacy_display()}</span>' )
    return label

@register.simple_tag
def job_color(job):
    """
    Returns a color based on job status.
    """
    return JOB_COLORS.get(job.state, "")


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
def echo(obj):
    # Makes debugging templates a bit easier
    print()

    print (dir(obj), "DEBUGGING")
    return ''


@register.filter
def can_edit(user, instance):
    """Returns true is instance is editable by user."""

    if user.is_superuser or instance.owner == user:
        return True

    return False

@register.filter
def can_create(user):
    """Returns true if user may create a new object"""

    return user.is_authenticated()

@register.inclusion_tag('widgets/project_header.html')
def project_header(project, label='Label', summary=True):
    """
    Returns a label for data sizes.
    """
    return dict(project=project, label=label, summary=summary)


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
    color = "" if data.data_type == const.GENERIC_TYPE else "green"
    label = const.DATA_TYPES.get(data.data_type, "Generic")
    return dict(label=label, color=color)


@register.inclusion_tag('widgets/form_nonfield_errors.html')
def form_nonfield_errors(form):
    """
    Turns the error lists into a dictionary that can be iterated over.
    """
    errorlist = list(form.non_field_errors())
    context = dict(errorlist=errorlist)
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


@register.filter(name='is_checkbox')
def is_checkbox(field):
    """
    Returns True if a field is a checkbox.
    """
    cond = isinstance(field, forms.BooleanField)
    return cond


@register.filter(name='is_selection')
def is_selection(field):
    """
    Returns True if a field's widget is a Selection
    """

    cond = isinstance(field.widget, widgets.Select)

    return cond
