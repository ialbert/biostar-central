import hjson as json
from pipeline.const import *
import uuid



def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

def join(*args):
    return os.path.abspath(os.path.join(*args))


def rewrite_specs(new_specs, file):

    return


def make_tmp_jsonfile(json_text, analysis_id):

    return


def check_fields(json_obj):

    required_keys = ["value", "display_type"]
    required_fields = ["template"]

    for field in required_fields:
        if field not in json_obj:
            raise Exception(f"'{field}' field is required in the json spec file.")
    for check in json_obj:
        data = json_obj[check]

        # Check Required keys
        for key in required_keys:
            if data.get(key) == None:
                raise KeyError(f"missing required key '{key}' in input: {data}")

        if data.get("label") == None and data.get("visible") == 1:
            data["label"] = check[0].upper() + check[1:]

    return json_obj


def safe_load(json_file):

    json_file = check_fields(json.load(open(json_file)))

    return json_file


def safe_loads(json_string):

    json_string = check_fields(json.loads(json_string))

    return json_string