from tempfile import TemporaryFile

import hjson
import logging
import uuid
from django.core.files import File

from . import tasks
from .const import *
from .models import Data, Analysis, Job, Project

CHUNK = 100

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


def make_toc(path):
    """
    Generate a table of contents into a temporary file.
    """

    size = 0

    def crawl(location, collect):
        nonlocal size
        for item in os.scandir(location):
            if item.is_dir():
                crawl(item.path, collect=collect)
            else:
                size += item.stat().st_size
                collect.append(os.path.abspath(item.path))

    lines = []
    crawl(path, collect=lines)
    text = '\n'.join(lines)
    fp = TemporaryFile()
    fp.write(text.encode('utf8'))
    fp.seek(0)
    return fp, lines, size


def create_data(project, user=None, stream=None, path=None, name=None, text='', summary='', data_type=None,
                link=False):
    size = 0

    # If the path is a directory, create the table of contents.
    if os.path.isdir(path):
        fp, lines, size = make_toc(path)
        link = False
        stream = File(fp)
        logger.info(f"Processing a directory.")
        name = f"Directory: {os.path.basename(fname)}"
        summary = f'Contains {len(lines)} files.'

    # The path is a file.
    if os.path.isfile(path):
        size = os.stat(path).st_size
        stream = File(open(path, 'rb'))
        name = os.path.basename(path)
        logger.info(f"Processing a file.")

    # Create the data.
    owner = user or project.owner
    name = name or "data.bin"
    data = Data.objects.create(name=name, owner=owner, state=Data.READY, text=text, project=project,
                               data_type=data_type, summary=summary)

    # We will allow invalid data to be added
    # while we build the site. TODO: be more strict here.
    if not stream:
        data.state = Data.PENDING
        data.save()
        link = True
        logger.error("Invalid stream specified.")
        # raise Exception(f"Empty stream. fname={path}")

    # Linking only points to an existing path
    if link:
        data.link = path
        data.save()
        logger.info(f"Linking to: {data.get_path()}")
    else:
        data.file.save(name, stream, save=True)
        logger.info(f"Saving to: {data.get_path()}")

    if data.can_unpack():
        logger.info(f"uwsgi active: {tasks.HAS_UWSGI}")
        if tasks.HAS_UWSGI:
            data_id = tasks.int_to_bytes(data.id)
            tasks.unpack(data_id=data_id).spool()
        else:
            tasks.unpack(data_id=data.id)

    # Updating data size.
    Data.objects.filter(pk=data.pk).update(size=size)

    logger.info(f"Added data id={data.pk}, type={data.data_type} name={data.name}")

    return data
