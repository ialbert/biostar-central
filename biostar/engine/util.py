import gzip, zipfile
import mimetypes
import os
import quopri
import tarfile
import uuid
from itertools import islice


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))

def fix_endings(text):
    return text.replace("\r\n", "\n")


class File(object):
    pass

def project_files(project):
    """
    Lists project data as they were files
    """

    data = project.data_set.order_by("sticky", "-date").all()

    def transform(f):
        b = File()
        b.path = f'{f.name}'
        b.is_dir = True
        b.name, b.size = f.name, f.size
        return b

    files = map(transform(data))

    return files

def scan_files(relpath, abspath, exclude=[]):
    """
    Generates a list of file objects at an absolute path.s
    """

    # Pathlike objects with attributes such as name, is_file
    files = list(filter(lambda p: p.name not in exclude, os.scandir(abspath)))

    # Sort the file list. Directories first, then by name.
    files = sorted(files, key=lambda p: (p.is_file(), p.name))

    def transform(f):
        b = File()
        b.path = f'{relpath}/{f.name}' if relpath else f.name
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
            #TODO:Should we keep doing this to files we do not know?
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


