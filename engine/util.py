import gzip
import mimetypes
import os
import quopri
import tarfile
import uuid
from itertools import islice

import hjson as json


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def smart_preview(fname):
    CHUNK_SIZE, LINE_COUNT = 1024, 10
    try:
        mimetype, mimecode = mimetypes.guess_type(fname)

        # print (mimetype, mimecode)

        if mimetype == 'application/x-tar' and mimecode == 'gzip':
            # A Tar gzip file
            tar = tarfile.open(fname, mode='r')
            lines = [f'{t.name}' for t in tar]
            text = "\n".join(lines)
        elif mimetype == None and mimecode == 'gzip':
            # A GZIP file.
            data = gzip.GzipFile(fname, 'r').read(CHUNK_SIZE)
            text = quopri.encodestring(data)
        elif mimetype == 'text/plain':
            stream = open(fname, 'rt')
            stream = islice(stream, LINE_COUNT)
            lines = [ line for line in stream]
            text = '\n'.join(lines)
        else:
            stream = open(fname, 'rb')
            data = stream.read(CHUNK_SIZE)
            text = quopri.encodestring(data)

    except Exception as exc:
        text = f'Preview error: {exc}'

    return text