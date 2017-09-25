from django import template

register = template.Library()


@register.filter(name='upper')
def upper(value):
    return value.upper()

@register.filter(name='date_sort')
def date_sort(instance):
    return instance.order_by("-id")

@register.inclusion_tag('includes/breadcrumb.html')
def breadcrumb(steps):
    return dict(steps=steps)


@register.inclusion_tag('includes/menubar.html')
def menubar(can_create_project=False, can_create_data=False,
            can_create_analysis=False):

    return dict(can_create_project=can_create_project,
                can_create_data=can_create_data,
                can_create_analysis=can_create_analysis)

