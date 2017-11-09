from django import forms
from .models import Data
from biostar.engine import const



def float_field(data):

    numrange = data.get("range", [1.0, 1000.0])
    min_value, max_value = numrange[0], numrange[1]

    label = data.get("label")
    widget = forms.NumberInput()
    help_text = data.get("help", f"Enter number between {min_value} and {max_value}")
    initial = data.get("value", 1)

    field = forms.FloatField(widget=widget, initial=initial, min_value=min_value, max_value=max_value,
                             help_text=help_text, label=label, required=False)

    return field



def select_field(data, choicefunc=None):

    if choicefunc:
        choices = choicefunc() or []
    else:
        choices = data.get("choices", [])

    initial = data.get("value", "")
    label = data.get("label", "")
    help_text = data.get("help", "")

    widget = forms.Select(choices=choices, attrs={"class":"ui dropdown"})
    field = forms.CharField(widget=widget, initial=initial, label=label, help_text=help_text)

    return field


def char_field(data):

    initial = data.get("value", "")
    label = data.get("label", "")
    help_text = data.get("help", "")

    field = forms.CharField(initial=initial, label=label, help_text=help_text)

    return field


def radioselect_field(obj):

    choices = obj.get("choices", [])
    initial = obj.get("value", "")
    label = obj.get("label", "")
    help_text = obj.get("help", "")

    widget = forms.RadioSelect(choices=choices)
    field = forms.CharField(widget=widget, initial=initial, label=label, help_text=help_text)

    return field


def number_field(data):

    numrange = data.get("range", [0, 1])
    min_value, max_value = min(numrange), max(numrange)
    label = data.get("label", "")
    widget = forms.NumberInput()
    help_text = data.get("help", f"Range: {min_value} and {max_value}")
    initial = data.get("value", 0)

    field = forms.IntegerField(
        label=label, initial=initial, min_value=min_value, max_value=max_value,
        help_text=help_text, widget=widget
    )

    return field


def file_field(data):
    widget = forms.FileInput()
    label = data.get("label", "")
    initial = data.get("value", "")
    field = forms.FileField(widget=widget, label=label, required=False, initial=initial)
    return field


def checkbox_field(data):
    label = data.get("label", "")
    help_text = data.get("help", "")
    initial = data.get("value", False)
    widget = forms.CheckboxInput

    field = forms.BooleanField(initial=initial, widget=widget, label=label, help_text=help_text, required=False)

    return field


def ignore(data):
    return ''


def data_field_generator(field, project, data_type=None):

    valid_type = const.DATA_TYPE_SYMBOLS.get(data_type, -1)

    query = Data.objects.filter(project=project, data_type=valid_type)

    datamap = dict((obj.id, obj) for obj in query)

    def choice_func():
        choices = [(d.id, d.name) for d in datamap.values()]
        return choices

    return select_field(field, choicefunc=choice_func)


TYPE2FUNC = {
    const.RADIO: radioselect_field,
    const.DROPDOWN: select_field,
    const.INTEGER: number_field,
    const.TEXTBOX: char_field,
    const.FLOAT: float_field,
    const.UPLOAD: file_field,
    const.CHECKBOX: checkbox_field,
}


