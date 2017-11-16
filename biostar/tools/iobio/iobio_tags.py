from django import template
from django.utils.safestring import mark_safe

register = template.Library()


@register.inclusion_tag('iobio/iobio.html')
def iobio(path, name):
    """
    Generates iobio bam.
    """
    return dict(path=path, name=name)
