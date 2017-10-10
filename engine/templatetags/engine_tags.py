from django import template
from django import forms
from engine import models

register = template.Library()


@register.filter(name='upper')
def upper(value):
    return value.upper()

@register.filter(name='date_sort')
def date_sort(instance):
    return instance.order_by("-id")


@register.filter(name='state_color')
def state_color(instance):

    states = {
              models.Job.QUEUED:{"color":"black",
                 "state": "Queued"},
              2:{"color":"yellow",
                 "state": "Running"},
              3:{"color":"green",
                 "state":"Finished"},
              4:{"color":"red",
                 "state":"Stopped"},
              }

    return states[instance]["color"]

# get rid of oth
@register.filter(name='int2state')
def int2state(instance):

    states = {
              1:"Queued",
              2:"Running",
              3:"Finished",
              4:"Stopped"
              }

    return states[instance]


@register.inclusion_tag('includes/breadcrumb.html')
def breadcrumb(steps):
    return dict(steps=steps)


@register.inclusion_tag('includes/menubar.html')
def menubar(can_create_project=False, can_create_data=False,
            can_create_analysis=False):

    return dict(can_create_project=can_create_project,
                can_create_data=can_create_data,
                can_create_analysis=can_create_analysis)

@register.filter(name='is_checkbox')
def is_checkbox(value):
    return isinstance(value, forms.BooleanField)