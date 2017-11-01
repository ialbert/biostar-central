import hjson
import logging
from django.core.files import File

from . import tasks
from .const import *
from .models import Data, Analysis, Job

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


def create_project(user, project_model):
    # project = project_model.objects.filter

    logger.info(f"{user.email} created project {project_model.id}")
    pass


def get_analysis(job=None,id=None, name=None, type=Analysis.USER):

    pass


def create_analysis(project, json_text, template, user=None, summary='', name='', text='', type=None):
    owner = user or project.owner
    name = name or 'Analysis name'
    text = text or 'Analysis text'
    type = type or Analysis.USER

    analysis = Analysis.objects.create(project=project, summary=summary, json_text=json_text,
                                       owner=owner, name=name, text=text, type=type,
                                       template=template)

    logger.info(f"Created analysis id={analysis.id} of type={analysis.get_type_display()}")

    return analysis


def edit_analysis():
    return


def create_job(analysis, user=None, project=None, json_text='', json_data={}, name=None, state=None, type=None):

    name = name or analysis.name
    state = state or Job.QUEUED
    owner = user or analysis.project.owner
    type = type or analysis.type
    project = project or analysis.project

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or analysis.json_text

    job = Job.objects.create(name=name, summary=analysis.summary, state=state, json_text=json_text,
                                   project=project, analysis=analysis, owner=owner, type=type,
                                   template=analysis.template)

    logger.info(f"Created job: '{job.name}'")
    logger.info(f"uwsgi active: {tasks.HAS_UWSGI}")
    if tasks.HAS_UWSGI:
        # Execute admin jobs asap
        if type == Job.ADMIN:
            jobid = (job.id).to_bytes(5, byteorder='big')
            tasks.execute_admin_job.spool(job_id=jobid)

    return job


def copy_data():
    return


def create_data(project, user=None, stream=None, fname=None, name="data.bin", text='', data_type=None, type=None):
    if fname:
        stream = File(open(fname, 'rb'))
        name = os.path.basename(fname)

    owner = user or project.owner
    text = text or "No description"
    data_type = data_type or GENERIC_TYPE
    type = type or Data.USER
    data = Data(name=name, owner=owner, type=type,
                text=text, project=project, data_type=data_type)

    # Need to save before uid gets triggered.
    data.save()
    # This saves the into the
    data.file.save(name, stream, save=True)

    if data.can_unpack():
        data.set_ready_flag(Data.PENDING)

        # Gets 1st admin analysis with unpack in its name and creates a job for that
        # can also take id.
        analysis = get_analysis(name="unpack", type=Analysis.ADMIN)

        create_job(analysis=analysis)

    else:
        data.set_ready_flag(Data.READY)

    # Updates its own size.
    data.set_size()

    logger.info(f"Added data id={data.id} of type={data.get_type_display()}")

    return data