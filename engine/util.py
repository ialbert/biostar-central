import hjson as json
import os
from .settings import BASE_DIR
from .factory import *
from pipeline.const import *


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


def join(*args):
    return os.path.abspath(os.path.join(*args))


def fill_makefile(filled_json, makefile_template):

    # will return a string.
    return NotImplemented


def rewrite_jsonspecs(new_specs, file):

    return


def remove_tmp_files():
    return


def make_tmp_jsonfile(json_text, analysis_id):


    tmp_jsonfile = join(BASE_DIR, '..', 'media',
                        f"tmp{analysis_id}.hjson")

    json_obj = json.loads(json_text)

    # write to tmp file
    json.dump(json_obj, open(tmp_jsonfile, "w"))

    return tmp_jsonfile



def safe_load(json_file):

    required_keys = ["name", "value" , "display_type"]

    json_file = json.load(open(json_file))

    for check in json_file:

        # Check Required keys
        for key in required_keys:

            if check.get(key)== None :
                raise KeyError(f"missing required key '{key}' in input: {check}")

        if not (check["display_type"] in TYPE2FUNC):
            print(check["display_type"])
            raise Exception(f"Incorrect value for display_type. Options are:{TYPE2FUNC.keys()}")

        if check.get("label")==None:
            check["label"] = check["name"][0].upper() + check["name"][1:]

    return json_file
