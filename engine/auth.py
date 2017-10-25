from biostar.tools import defaults
from django.core.files import File
import hjson
import logging
from .const import *
import os

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

    #project = project_model.objects.filter

    logger.info(f"{user.email} created project {project_model.id}")
    pass


def create_analysis(user, project, analysis_model, json_text, template, summary='', name='', text='', type=None):

    owner = user or project.owner
    name = name or 'Analysis name'
    text = text or 'Analysis text'
    type = type or defaults.USAGE

    analysis = analysis_model.objects.create(project=project, summary=summary, json_text=json_text,
                                       owner=owner, name=name, text=text, type=type,
                                       template=template)
    return analysis


def edit_analysis():
    return



def create_data(user, data_model, project, stream=None, fname=None, name="data.bin", text='', data_type=None, type=None):

    if fname:
        stream = File(open(fname, 'rb'))
        name = os.path.basename(fname)

    owner = user or project.owner
    text = text or "No description"
    data_type = data_type or GENERIC_TYPE
    type = type or defaults.USAGE
    data = data_model(name=name, owner=owner, type=type,
                text=text, project=project, data_type=data_type)

    # Need to save before uid gets triggered.
    data.save()
    # This saves the into the
    data.file.save(name, stream, save=True)

    # Set the pending to ready after the file saves.
    data.ready_state()

    # Updates its own size.
    data.set_size()
    return data

# Can we put stuff from views in here?
def edit_data():
    return


def create_job(user, analysis, job_model,project=None, json_text='', json_data={}, name=None, state=None, type=None):

    name = name or analysis.name
    state = state or job_model.QUEUED
    owner = user or analysis.project.owner
    type = type or defaults.USAGE
    project = project or analysis.project

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or analysis.json_text

    job = job_model.objects.create(name=name, summary=analysis.summary, state=state, json_text=json_text,
                             project=project, analysis=analysis, owner=owner, type=type,
                             template=analysis.template)
    return job

