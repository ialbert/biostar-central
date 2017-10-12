from django import forms
from .models import Data
from pipeline.const import *


def float_field(field):

    numrange = field.get("range", [1.0,1000.0])
    min_value, max_value = numrange[0], numrange[1]

    label = field.get("label")
    widget = forms.NumberInput()
    help_text = field.get("help", f"Enter number between {min_value} and {max_value}")
    initial = field.get("value", 1)

    field = forms.FloatField(widget=widget,
                               initial=initial,
                               min_value=min_value,
                               max_value=max_value,
                               help_text=help_text,
                               label=label,
                               required=False)

    return field


def select_field(field):

    choices = get_choices(field)
    initial = field.get("value")
    label = field.get("label")
    help_text = field.get("help", "Pick values from a dropdown menu")

    widget = forms.Select(choices=choices)

    field = forms.CharField(widget=widget,
                            initial=initial,
                            label=label,
                            help_text=help_text)

    return field


def radioselect_field(field):

    choices = get_choices(field)

    initial = field.get("value")
    label = field.get("label")
    help_text = field.get("help", "Enter one of more option to display")

    widget = forms.RadioSelect(choices=choices)

    field = forms.CharField(widget=widget,
                            initial=initial,
                            label=label,
                            help_text=help_text)

    return field


def number_field(field):

    numrange = field.get("range", [1,1000])

    min_value, max_value = numrange[0], numrange[1]
    label = field.get("label")
    widget = forms.NumberInput()
    help_text = field.get("help", f"Enter number between {min_value} and {max_value}")
    initial = field.get("value", 1)

    field = forms.IntegerField(widget=widget,
                               initial=initial,
                               min_value=min_value,
                               max_value=max_value,
                               help_text=help_text,
                               label=label)

    return field


def file_field(field):

    widget = forms.FileInput()
    label = field.get("label")
    initial = field.get("value")

    field = forms.FileField(widget=widget,
                            label=label,
                            required=False,
                            initial=initial)
    return field


def checkbox_field(field):

    boolmap = {'true':True, 'false':False}

    label = field.get("label")
    help_text = field.get("help", "Check option for true.")

    inital = boolmap.get(field.get("value", "false"))
    widget = forms.CheckboxInput

    field = forms.BooleanField(initial=inital,
                               widget=widget,
                               label=label,
                               help_text=help_text,
                               required=False)
    return field


def model_specs(field):
    return


def handle_scripts(field):
    return


def get_choices(field):

    # View data already in database
    choices = field.get("choices")
    if choices:
        choices = list(choices.items())

    if field.get("origin") == "PROJECT":
        # Just loads all data for now ( not project specific ).
        data = Data.objects.all()
        choices = []
        for d in data:
            choices.append((d.id, d.title))

    return choices


TYPE2FUNC = {

    RADIO: radioselect_field,
    DROPDOWN: select_field,
    INTEGER: number_field,
    FLOAT : float_field,
    UPLOAD : file_field,
    CHECKBOX :checkbox_field,
    MODEL: model_specs,
    SCRIPT : handle_scripts,

}

def check_display(display_type):

    if not (display_type in TYPE2FUNC):
        raise Exception(f"Incorrect value for display_type. Options are:{TYPE2FUNC.keys()}")

    return
