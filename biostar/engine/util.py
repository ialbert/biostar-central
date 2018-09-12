import gzip
import io
import mimetypes
import os
import quopri
import tarfile
import uuid
from itertools import islice
from urllib.parse import quote

CHUNK = 1024 * 1024


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def fix_endings(text):
    return text.replace("\r\n", "\n")


class File(object):
    pass


FILE_ICON = '''
<div class="item">
  <i class="file icon"></i> %s
</div>
'''

FOLDER_ICON = '''
<div class="item">
  <i class="blue folder icon"></i> %s
  <div class="list">
'''

FILE_ICON = FILE_ICON.strip()
FOLDER_ICON = FOLDER_ICON.strip()

def directory_tree(path, collect=[]):
    for entry in os.scandir(path):
        if entry.is_dir():
            collect.append(FOLDER_ICON % entry.name)
            directory_tree(entry, collect=collect)
            collect.append(f'</div>')
            collect.append(f'</div>')
        else:
            collect.append(FILE_ICON % entry.name)
    return collect

def scan_files(relpath, abspath, root, exclude=[]):
    """
    Generates a list of file objects at an absolute path.s
    """
    if not (abspath.startswith(root) and os.path.exists(abspath)):
        raise Exception("Can not access an invalid directory.")

    # Pathlike objects with attributes such as name, is_file
    files = list(filter(lambda p: p.name not in exclude, os.scandir(abspath)))

    # Sort the file list. Directories first, then by name.
    files = sorted(files, key=lambda p: (p.is_file(), p.name))

    def transform(f):
        b = File()
        b.path = os.path.join(relpath, f.name) if relpath else f.name
        b.is_dir = f.is_dir()
        if os.path.exists(f.path):
            b.name, b.size = f.name, f.stat().st_size
        else:
            b.name, b.size = f.name, 0
        return b

    files = map(transform, files)
    files = list(files)

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
    mode = 'w' if isinstance(stream, io.StringIO) else 'wb'

    with open(dest, mode) as fp:
        chunk = stream.read(CHUNK)
        while chunk:
            fp.write(chunk)
            chunk = stream.read(CHUNK)

    return dest


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
