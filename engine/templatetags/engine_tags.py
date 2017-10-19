from django import forms
from django import template
from engine.models import Job, make_html
from django.utils.safestring import mark_safe
from textwrap import dedent
register = template.Library()

JOB_SEGMENT_COLORS = {
    Job.ERROR: "red", Job.QUEUED: "blue", Job.RUNNING: "orange", Job.FINISHED:"green"
}


@register.simple_tag
def job_color(job):
    return JOB_SEGMENT_COLORS.get(job.state, "")


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


@register.filter(name='date_sort')
def date_sort(instance):
    return instance.order_by("-id")


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


@register.filter(name='convert_bytes')
def convert_bytes_to(bytes, to="mega"):
    a = {'kilo': 1, 'mega': 2, 'giga': 3, 'tera': 4, 'penta': 5}
    r = float(bytes)
    for i in range(a[to]):
        r = r / 1000

    return f"{r:.2f}"
