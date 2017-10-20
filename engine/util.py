import hjson as json
import os, mimetypes, quopri, gzip, tarfile

from biostar.tools.const import *
import uuid


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]

def join(*args):
    return os.path.abspath(os.path.join(*args))

def smart_preview(fname):
    CHUNK_SIZE, LINE_COUNT = 1024, 10
    try:
        mimetype, mimecode = mimetypes.guess_type(fname)

        #print (mimetype, mimecode)

        stream = open(fname, 'rb')

        if mimetype == 'application/x-tar' and mimecode == 'gzip':
            # A Tar gzip file
            tar = tarfile.open(fname, mode='r')
            lines = [ f'{t.name}' for t in tar ]
            text = "\n".join(lines)
        elif mimetype == None and mimecode == 'gzip':
            # A GZIP file.
            data = gzip.GzipFile(fname, 'r').read(CHUNK_SIZE)
            text = quopri.encodestring(data)
        elif mimetype == 'text/plain':
            lines = [next(stream) for x in range(LINE_COUNT)]
            text = '\n'.join(lines)
        else:
            data = stream.read(CHUNK_SIZE)
            text = quopri.encodestring(data)

    except Exception as exc:
        text = f'Preview error: {exc}'

    return text

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
