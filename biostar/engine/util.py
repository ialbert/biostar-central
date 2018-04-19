import gzip, zipfile
import io
import mimetypes
import os
import quopri
import tarfile
import uuid
from itertools import islice


CHUNK = 1024 * 1024

class InvalidDirectoryError(Exception):
    def __str__(self):
        return "Can not access an invalid directory."


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))

def fix_endings(text):
    return text.replace("\r\n", "\n")


class File(object):
    pass


def scan_files(relpath, abspath, root, exclude=[]):
    """
    Generates a list of file objects at an absolute path.s
    """
    if not (abspath.startswith(root) and os.path.exists(abspath)) :
        raise InvalidDirectoryError

    # Pathlike objects with attributes such as name, is_file
    files = list(filter(lambda p: p.name not in exclude, os.scandir(abspath)))

    # Sort the file list. Directories first, then by name.
    files = sorted(files, key=lambda p: (p.is_file(), p.name))

    def transform(f):
        b = File()
        b.path = os.path.join(relpath, f.name) if relpath else f.name
        b.is_dir = f.is_dir()

        b.name, b.size = f.name, f.stat().st_size
        return b

    files = map(transform, files)

    return files

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
            lines = [ line.strip() for line in stream]
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

    mode = 'w' if isinstance(stream , io.StringIO) else 'wb'

    with open(dest, mode) as fp:
        chunk = stream.read(CHUNK)
        while chunk:
            fp.write(chunk)
            chunk = stream.read(CHUNK)

    return dest



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
