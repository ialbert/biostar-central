from django import forms
from django import template
from biostar.engine.models import Job, make_html
from django.utils.safestring import mark_safe
from textwrap import dedent
from biostar.engine import const
from django.contrib.staticfiles.templatetags.staticfiles import static

register = template.Library()

JOB_SEGMENT_COLORS = {
    Job.ERROR: "red", Job.QUEUED: "blue", Job.RUNNING: "teal", Job.FINISHED: "green"
}


@register.simple_tag
def job_color(job):
    return JOB_SEGMENT_COLORS.get(job.state, "")

@register.simple_tag
def img(obj):

    if obj.image:
        return obj.image.url
    else:
        return static("images/placeholder.png")


@register.simple_tag
def type_label(data):
    label = const.DATA_TYPES.get(data.data_type, "Generic")
    return label

@register.inclusion_tag('widgets/form_nonfield_errors.html')
def form_nonfield_errors(form):
    errorlist = list(form.non_field_errors())
    context = dict(errorlist=errorlist)
    return context

@register.simple_tag
def field_state(field):
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
def is_checkbox(value):
    return isinstance(value, forms.BooleanField)

