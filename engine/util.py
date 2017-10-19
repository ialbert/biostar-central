import hjson as json
import os
from biostar.tools.const import *
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

    return json_obj


def safe_load(json_file):

    json_file = check_fields(json.load(open(json_file)))

    return json_file


def safe_loads(json_string):

    json_string = check_fields(json.loads(json_string))

    return json_string
