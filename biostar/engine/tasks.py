'''
This is active only when deployed via UWSGI
'''

import logging, time, shutil, subprocess
from django.core import management
from django.utils.text import force_text

logger = logging.getLogger("engine")

HAS_UWSGI = False

def int_to_bytes(value):
    return (value).to_bytes(8, byteorder='big')

def int_from_bytes(args, name):
    return int.from_bytes(args[name].encode(), byteorder='big')

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @timer(3)
    def execute_timer(args):
        from biostar.engine.models import Job
        # logger.info(f"executing timer with {args}")
        # It is faster to check here
        first = Job.objects.filter(state=Job.QUEUED).first()
        if first:
            management.call_command('job', next=True)

    @spool
    def execute_job(args):
        '''
        Spools at normal priority.
        '''
        job_id = int_from_bytes(args, "job_id")
        management.call_command('job', id=job_id)

    @spool
    def unpack(args):
        data_id = int_from_bytes(args, "data_id")
        unpacker(data_id=data_id)

    @spool
    def copy(args):
        return

except ModuleNotFoundError as exc:

    # No spooling same interface.
    def unpack(data_id):
        unpacker(data_id=data_id)


def execute(command, workdir="."):
    proc = subprocess.run(command, cwd=workdir, shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = force_text(proc.stdout), force_text(proc.stderr)
    return stdout, stderr


def unpacker(data_id):
    """
    Unpacks a data and sets the status to READY.
    """
    from biostar.engine.models import Data
    query = Data.objects.filter(id=data_id)
    data = query.first()
    query.update(state=Data.PENDING)
    command = f'tar xzvf {data.file.path}'
    try:
        stdout, stderr = execute(command=command, workdir=data.get_datadir())
        logger.info(f'Data id={data_id} has been unpacked')
        query.update(state=Data.READY)
    except Exception as exc:
        logger.error(f"Error: f{exc}")
        query.update(state=Data.ERROR)


def copier(source=None, target_data=None, target_project=None, fname=None):

    """
    Copies source data to target_data id. or adds a fname to target_project id
    """
    from biostar.engine.models import Data, Project
    from biostar.engine import auth

    assert (source and target_data) or (target_project and fname)

    #TODO: should probs link em instead
    source = Data.objects.filter(id=source).first()
    target_data = Data.objects.filter(id=target_data).first()

    project = Project.objects.filter(id=target_project).first()

    if project:
        assert fname
        auth.create_data(project=project, fname=fname)
        return
    if source:
        assert target_data

    #TODO:Test this part
    shutil.copy(source.file.path, target_data.file.path)
    return