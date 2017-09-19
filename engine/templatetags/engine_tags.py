from django import template

register = template.Library()


@register.filter(name='upper')
def upper(value):
    return value.upper()

@register.simple_tag
def foo(value):
    return f"FOO: {value}"


@register.inclusion_tag('includes/breadcrumb.html')
def breadcrumb(steps):
    return dict(steps=steps)

@register.inclusion_tag('includes/menubar.html')
def menubar(action):
    return dict(action=action)
