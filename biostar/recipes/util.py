import gzip
import io
import mimetypes
import os
import quopri
import tarfile
import uuid
import bleach
import shlex
import random
from itertools import islice
from urllib.parse import quote
from datetime import datetime
from django.utils.timezone import utc
import toml as hjson

CHUNK = 1024 * 1024


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def fix_endings(text):
    return text.replace("\r\n", "\n")


class File(object):
    pass


def now():
    return datetime.utcnow().replace(tzinfo=utc)



def pp(json):
    """
    Pretty prints json
    """
    text = hjson.dumps(json)
    return text


def smart_preview(fname):
    CHUNK_SIZE, LINE_COUNT = 1024, 10
    try:
        mimetype, mimecode = mimetypes.guess_type(fname)

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
            lines = [line.strip() for line in stream]
            text = '\n'.join(lines)
        else:
            try:
                stream = open(fname, 'rt')
                text = stream.read(CHUNK_SIZE)
            except Exception as exc:
                stream = open(fname, 'rb')
                data = stream.read(CHUNK_SIZE)
                text = quopri.encodestring(data)

    except Exception as exc:
        text = f'Preview error: {exc}'

    return text


def write_stream(stream, dest):

    # Output needs to be opened in "w" if the incoming stream is a string
    # and "wb" if the incoming stream is a file.
    mode = "w" if isinstance(stream, io.StringIO) else "wb"

    with open(dest, mode) as fp:
        chunk = stream.read(CHUNK)
        while chunk:
            fp.write(chunk)
            chunk = stream.read(CHUNK)
    return dest


def clean_text(value):
    #TODO: investigate more,
    # shoule be applied only to bash scripts.
    return value
    #return shlex.quote(textbox)


def qiime2view_link(file_url):
    template = "https://view.qiime2.org/visualization/?type=html&src="

    file_url = quote(string=file_url, safe="")

    return template + file_url



def findfiles(location, collect):
    """
    Returns a list of all files in a directory.
    """

    for item in os.scandir(location):

        if item.is_dir():
            findfiles(item.path, collect=collect)
        else:
            collect.append(os.path.abspath(item.path))

    return collect
