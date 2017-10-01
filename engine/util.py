import json
from .settings import ALLOWED_WIDGETS


#FORM_TYPES = []


def fill_makefile(filled_json, makefile_template):

    # will return a string.
    return NotImplemented



def safe_load(json_file):

    required_keys = ["name", "value" ]

    json_file = json.load(json_file)

    for check in json_file:

        # Check Required keys
        for key in required_keys:

            if check.get(key)== None :
                raise KeyError(f"{check} missing required key {key}")

        # Verify widget
        if check.get("widget"):
            assert check["widget"] in ALLOWED_WIDGETS
            # Get correct from type from widget
            check["form_type"] = ALLOWED_WIDGETS[check["widget"]]

            if check["form_type"] == "IntegerField":
                assert check.get("range")
                check["min_value"], check["max_value"] = check["range"][0], check["range"][1]

        if check.get("label")==None:
            check["label"] = check["name"][0].upper() + check["name"][1:]

    return json_file
