from django import template

register = template.Library()


@register.filter(name='upper')
def upper(value):
    return value.upper()

@register.inclusion_tag('includes/breadcrumb.html')
def breadcrumb(steps):
    return dict(steps=steps)


@register.inclusion_tag('includes/menubar.html')
def menubar(is_project_list=False, is_data_list=False,
            is_analysis_list=False):

    return dict(is_project_list=is_project_list,
                is_data_list=is_data_list,
                is_analysis_list=is_analysis_list)

