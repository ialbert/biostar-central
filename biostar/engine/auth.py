import hjson, logging, shutil, tarfile
from django.core.files import File
from . import tasks
from .const import *
from .models import Data, Analysis, Job, Project


logger = logging.getLogger("engine")


def get_data(user, project, query, data_type=None):
    """
    Returns a dictionary keyed by data stored in the project.
    """
    if data_type:
        query = query.filter(data_type=data_type)
    datamap = dict((obj.id, obj) for obj in query)
    logger.info(f"{user.email} got {len(datamap)} data from {project.name}")

    return datamap


def create_project(user, name, uid='', summary='', text=''):

    print (uid)
    project = Project.objects.create(
        name=name, uid=uid,  summary=summary, text=text, owner=user,
    )
    logger.info(f"Created project: {project.name} uid: {project.uid}")
    pass


def create_analysis(project, json_text, template,
                    uid=None, user=None, summary='', name='', text=''):
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

def create_data(project, user=None, stream=None, fname=None, name="data.bin", text='', data_type=None):

    if fname:
        stream = File(open(fname, 'rb'))
        name = os.path.basename(fname)

    owner = user or project.owner
    text = text or "No description"
    data_type = data_type or GENERIC_TYPE

    # Create the data
    data = Data.objects.create(name=name, owner=owner, state=Data.READY,
                               text=text, project=project, data_type=data_type)

    # This saves the into the
    data.file.save(name, stream, save=True)

    if data.can_unpack():

        logger.info(f"uwsgi active: {tasks.HAS_UWSGI}")
        if tasks.HAS_UWSGI:
            data_id = tasks.int_to_bytes(data.id)
            tasks.unpack(data_id=data_id).spool()
        else:
            tasks.unpack(data_id=data.id)

    # Updates its own size.
    data.set_size()

    logger.info(f"Added data id={data.name}, name={data.name}")

    return data
