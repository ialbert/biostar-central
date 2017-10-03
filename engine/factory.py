
from string import Template
from .models import Data


FULL_TEMPLATE = r"""self.fields['json_'+ '$name']= forms.$form_type(
                                             widget=forms.$widget(choices=$choices),
                                             initial='$value',
                                             label='$label',
                                             min_value=$min_value,
                                             max_value=$max_value )"""
def selectfield(field):

    global FULL_TEMPLATE

    choices = get_choices(field)

    data = {"name": field["name"], "form_type": field["form_type"],
            "label": field["label"], "value": field["value"],
            "widget": field["widget"], "choices":choices}

    template = FULL_TEMPLATE.replace("min_value=$min_value,", "")
    template = template.replace("max_value=$max_value", "")
    field_template = Template(template).safe_substitute(data)

    return field_template


def radiofield(field):

    global FULL_TEMPLATE

    choices = get_choices(field)

    data = {"name": field["name"], "form_type": field["form_type"],
            "label": field["label"], "value": field["value"],
            "widget": field["widget"], "choices": choices}

    template = FULL_TEMPLATE.replace("min_value=$min_value,", "")
    template = template.replace("max_value=$max_value", "")
    field_template = Template(template).safe_substitute(data)

    return field_template




def numberfield(field):

    global FULL_TEMPLATE

    data = {"name": field["name"], "form_type": field["form_type"],
            "label": field["label"], "value": field["value"],
            "widget": field["widget"], "max_value":field["max_value"],
            "min_value":field["min_value"]}

    template = FULL_TEMPLATE.replace("choices=$choices", "")
    field_template = Template(template).safe_substitute(data)

    return field_template


def filefield(field):

    global FULL_TEMPLATE

    data = {"name": field["name"], "form_type": field["form_type"],
            "label": field["label"], "value": field["value"],
            "widget": field["widget"]}

    template = FULL_TEMPLATE.replace("choices=$choices", "")
    template = template.replace("min_value=$min_value,", "")
    template = template.replace("max_value=$max_value", "")
    field_template = Template(template).safe_substitute(data)

    return field_template


def get_choices(feild):

    # View data already in database
    choices = feild.get("choices")

    if feild["name"] == "data" and feild.get("origin") == "PROJECT":
        # Just loads all data for now ( not project specific ).
        data = Data.objects.all()
        choices = []
        for d in data:
            choices.append((d.id, d.title))

    return choices