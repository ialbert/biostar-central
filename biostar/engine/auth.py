
import hjson
import logging
import uuid, shutil
from django.core.files import File
from django.utils.text import slugify

from . import tasks
from .const import *
from .models import Data, Analysis, Job, Project
from . import models

CHUNK = 1024 * 1024

logger = logging.getLogger("engine")

def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def get_project_list(user):
    """
    Return projects with privileges relative to a user.
    """
    query = Project.objects.all()

    # Superusers see everything
    if user.is_superuser:
        return query

    # Unauthenticated users see public projects.
    if user.is_anonymous:
        return query.filter(privacy=Project.PUBLIC)

    # We need to rework this first to allow adding users to groups.

    # get the private and sharable projects belonging to the same user
    # then merge that with the public projects query
    # query = query.filter(
    #              Q(owner=user),
    #              Q(privacy=Project.PRIVATE)|
    #              Q(privacy=Project.SHAREABLE)) | query.filter(privacy=Project.PUBLIC)

    return query


def get_data(user, project, query, data_type=None):
    """
    Returns a dictionary keyed by data stored in the project.
    """
    if data_type:
        query = query.filter(data_type=data_type)
    datamap = dict((obj.id, obj) for obj in query)
    logger.info(f"{user.email} got {len(datamap)} data from {project.name}")

    return datamap


def create_project(user, name, uid='', summary='', text='', stream='', privacy=Project.PRIVATE, sticky=True):
    project = Project.objects.create(
        name=name, uid=uid, summary=summary, text=text, owner=user, privacy=privacy, sticky=sticky)

    if stream:
        project.image.save(stream.name, stream, save=True)

    logger.info(f"Created project: {project.name} uid: {project.uid}")

    return project


def create_analysis(project, json_text, template, uid=None, user=None, summary='', name='', text=''):
    owner = user or project.owner
    name = name or 'Analysis name'
    text = text or 'Analysis text'

    analysis = Analysis.objects.create(project=project, uid=uid, summary=summary, json_text=json_text,
                                       owner=owner, name=name, text=text,
                                       template=template)

    logger.info(f"Created analysis: uid={analysis.uid}")

    return analysis


def create_job(analysis, user=None, project=None, json_text='', json_data={}, name=None, state=None, type=None):
    name = name or analysis.name
    state = state or Job.QUEUED
    owner = user or analysis.project.owner

    project = project or analysis.project

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or analysis.json_text

    job = Job.objects.create(name=name, summary=analysis.summary, state=state, json_text=json_text,
                             project=project, analysis=analysis, owner=owner,
                             template=analysis.template)

    logger.info(f"Created job: {job.name}")

    return job


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


def create_data_path(data, fname):

    # Limit the file name lenght.
    fname = os.path.basename(fname)[-100:]

    # Slugify each element of the path to make them "safe".
    pieces = map(slugify, fname.split("."))

    # Put it back together.
    fname = ".".join(pieces)

    # This will store the data.
    data_dir = data.get_data_dir()

    # Make the data directory if it does not exist.
    os.makedirs(data_dir, exist_ok=True)

    # Build file with original name.
    data_path = os.path.join(data_dir, fname)

    return data_path

def create_data(project, user=None, stream=None, path=None, name='data.bin',
                text='', summary='', data_type=None, link=False):

    # Create the data.
    data = Data.objects.create(name=name, owner=user, state=Data.PENDING,  project=project,
                               data_type=data_type, summary=summary, text=text)



    # If the path is a directory, symlink all files.
    if path and os.path.isdir(path):
        collect = findfiles(path, collect=[])
        for src in collect:
            dest = create_data_path(data, src)
            os.symlink(src, dest)
        logger.info(f"Symlinked {len(collect)} files.")

    # The path is a file and should be linked.
    if path and os.path.isfile(path):
        dest = create_data_path(data, path)
        if link:
            os.symlink(path, dest)
            logger.info(f"Symlinked to: {path}")
        else:
            logger.info(f"Copied file to: {dest}")
            shutil.copy(path, dest)

    missing = not (os.path.isdir(path) or os.path.isfile(path))
    if path and missing:
        logger.error(f"File not found for: {path}")

    # An incoming stream written into the destination.
    if stream:
        dest = create_data_path(data, path)
        with open(dest, 'wb') as fp:
            chunk = stream.read(CHUNK)
            while chunk:
                fp.write(chunk)
                chunk = stream.read(CHUNK)

    # Find all files in the data directory.
    collect = findfiles(data.get_data_dir(), collect=[])

    # Remove the table of contents from the collected files.
    collect.remove(data.get_path())

    # Write the table of contents.
    with open(data.file, 'wt') as fp:
        fp.write("\n".join(collect))

    # Finalize the data.
    Data.objects.filter(pk=data.pk).update(state=Data.READY, size=100)

    logger.info(f"Added data id={data.pk}, type={data.data_type} name={data.name}")

    return data
