from textwrap import dedent

from django import forms
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.forms import widgets
from django.utils.safestring import mark_safe

from biostar.engine import const
from biostar.engine.models import Job, make_html

register = template.Library()

JOB_COLORS = {
    Job.ERROR: "red", Job.QUEUED: "blue", Job.RUNNING: "teal", Job.COMPLETED: "green"
}

@register.simple_tag
def data_color(data):
    """
    Returns a color based on job status.
    """
    return "" if data.data_type == const.GENERIC_TYPE else "green"


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


@register.simple_tag
def type_label(data):
    """
    Returns readable names for data types.
    """
    label = const.DATA_TYPES.get(data.data_type, "Generic")
    return label


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
