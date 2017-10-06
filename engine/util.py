import hjson as json
#from .settings import BASE_DIR
from pipeline.const import *
import uuid



def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def fill_makefile(filled_json, makefile_template):

    # will return a string.
    return NotImplemented


def rewrite_jsonspecs(new_specs, file):

    return


def remove_tmp_files():

    return


def group_by_state(jobs_list):

    # return a dictionary with keys being
    # a state and value is a list of jobs in that state.
    return


def make_tmp_jsonfile(json_text, analysis_id):


    tmp_jsonfile = join(BASE_DIR, '..', 'media',
                        f"tmp{analysis_id}.hjson")

    json_obj = json.loads(json_text)

    # write to tmp file
    json.dump(json_obj, open(tmp_jsonfile, "w"))

    return tmp_jsonfile


def safe_load(json_file):

    required_keys = ["value" , "display_type"]

    json_file = json.load(open(json_file))

    for check in json_file:
        data = json_file[check]
        # Check Required keys
        for key in required_keys:

            if data.get(key)== None :
                raise KeyError(f"missing required key '{key}' in input: {data}")

        #if not (data["display_type"] in TYPE2FUNC):
        #    raise Exception(f"Incorrect value for display_type. Options are:{TYPE2FUNC.keys()}")

        if data.get("label")==None:
            data["label"] = check[0].upper() + check[1:]

    return json_file
