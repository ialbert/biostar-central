from django import template
from django.utils.safestring import mark_safe

register = template.Library()

@register.inclusion_tag('igv/bam.xml')
def bam(path, name):
    """
    Generates tracks for bam files.
    """
    return dict(path=path, name=name)

@register.inclusion_tag('igv/resources.xml')
def resources(path):
    """
    Generates a resource tag.
    """
    return dict(path=path)
