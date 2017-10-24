from biostar.tools import defaults
from django.core.files import File
import logging
from .const import *
import os

logger = logging.getLogger("engine")


def get_data(user, project, data_type, query):
    """
    Returns a dictionary keyed by data stored in the project.
    """
    if data_type:
        query = query.filter(data_type=data_type)
    datamap = dict((obj.id, obj) for obj in query)
    logger.info(f"{user.email} got {len(datamap)} data from {project.name}")

    return datamap


def create_project(user, project_model):

    #project = Project.objects.filter

    print("CREATED PROJECT TESTED")

    logger.info(f"{user.email} created project {project_model.id}")
    pass


def create_analysis(user, project, analysis_model, json_text, template, summary='', name='', text='', usage=None):

    owner = user or project.owner
    name = name or 'Analysis name'
    text = text or 'Analysis text'
    usage = usage or defaults.USAGE

    analysis = analysis_model.objects.create(project=project, summary=summary, json_text=json_text,
                                       owner=owner, name=name, text=text, usage=usage,
                                       template=template)
    return analysis


def edit_analysis():
    return



def create_data(user, data_model, project, stream=None, fname=None, name="data.bin", owner=None, text='', data_type=None, usage=None):

    if fname:
        stream = File(open(fname, 'rb'))
        name = os.path.basename(fname)

    owner = owner or project.owner
    text = text or "No description"
    data_type = data_type or GENERIC_TYPE
    usage = usage or defaults.USAGE
    data = data_model(name=name, owner=owner, usage=usage,
                text=text, project=project, data_type=data_type)

    # Need to save before uid gets triggered.
    data.save()
    # This saves the into the
    data.file.save(name, stream, save=True)

    # Set the pending to ready after the file saves.
    data_model.objects.filter(id=data.id).update(state=data_model.READY)

    # Updates its own size.
    data.set_size()
    return data



def edit_data():
    return



def create_job(user, data):

    name = name or self.name
    state = state or Job.QUEUED
    owner = owner or self.project.owner
    usage = usage or defaults.USAGE

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or self.json_text

    job = Job.objects.create(name=name, summary=self.summary, state=state, json_text=json_text,
                             project=self.project, analysis=self, owner=owner, usage=usage,
                             template=self.template)
    pass

