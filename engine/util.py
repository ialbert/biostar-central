import json
from .settings import ALLOWED_WIDGETS
from .factory import *


TYPE2FUNC = {

    "RADIO": radiofield,
    "SELECT": selectfield,
    "NUMBER": numberfield,
    "FILE" : filefield,
    "EMPTY":"",

}


def fill_makefile(filled_json, makefile_template):

    # will return a string.
    return NotImplemented


def handle_no_type(field):

    # "Invisibile" fields are handled more downstream
    if field.get("visible") != 1:
        return

    # do nothing when "type" is not a field key.
    return


def safe_load(json_file):

    required_keys = ["name", "value" , "type"]

    json_file = json.load(open(json_file))

    for check in json_file:

        # Check Required keys
        for key in required_keys:

            if check.get(key)== None :
                raise KeyError(f"{check} missing required key <'{key}'>")

        # Verify widget
        if check.get("widget"):
            if not (check["widget"] in ALLOWED_WIDGETS):
                raise Exception(f"{check['widget']} is not part of allowed widgets ({ALLOWED_WIDGETS.keys()})")

            # Map a form to the given widget
            check["form_type"] = ALLOWED_WIDGETS[check["widget"]]

            if check["form_type"] == "IntegerField":
                assert check.get("range")
                check["min_value"], check["max_value"] = check["range"][0], check["range"][1]

        if check.get("label")==None:
            check["label"] = check["name"][0].upper() + check["name"][1:]

    return json_file
