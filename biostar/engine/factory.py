from django import forms

from biostar.engine import const
from . import models

# Share the logger with models.
logger = models.logger













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

    widget = forms.Select(choices=choices, attrs={"class": "ui dropdown"})
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
    """
    Generates a SELECT field populated by data names that
    are of a certain type.
    """

    # Figure out the type from the symbol.
    valid_type = const.DATA_TYPE_SYMBOLS.get(data_type, -1)

    # Get the data that match this type in the project.
    query = models.Data.objects.filter(project=project, data_type=valid_type).order_by("sticky", "-date")

    # Create a mapping of data to id.
    datamap = dict((obj.id, obj) for obj in query)

    # The choice generator.
    def choice_func():
        choices = [(d.id, d.name) for d in datamap.values()]
        return choices

    # Returns a SELECT field with the choices.
    return select_field(field, choicefunc=choice_func)


def dynamic_field(data, project=None):
    """
    Creates a DJANGO form field from a dictionary.
    """

    # Get the known field types.
    field_types = get_field_types()

    if not hasattr(data, 'get'):
        # Not a "dictionary-like" data
        return None

    # Find out the display type.
    display_type = data.get("display_type")

    # Fields with no display type are not visible.
    if not display_type:
        return None

    # Data is accessed via paths or links.
    is_data = data.get("path") or data.get("link")

    if is_data and project:
        # Project specific data needs to be generated from known data.
        data_type = data.get("data_type")
        field = data_field_generator(data, project=project, data_type=data_type)
    else:
        # In all other cases we generate a field from the tupe.
        func = field_types.get(display_type)
        if not func:
            logger.error(f"Invalid display_type={display_type}")
            return None
        field = func(data)

    return field


def get_field_types():
    """
    Maps strings constants to field types.
    """
    field_types = {
        const.RADIO: radioselect_field,
        const.DROPDOWN: select_field,
        const.INTEGER: number_field,
        const.TEXTBOX: char_field,
        const.FLOAT: float_field,
        const.UPLOAD: file_field,
        const.CHECKBOX: checkbox_field,
    }

    return field_types
