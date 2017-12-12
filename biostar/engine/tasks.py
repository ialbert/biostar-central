'''
This is active only when deployed via UWSGI
'''

import logging, time, shutil, subprocess
from django.core import management
from django.utils.text import force_text
import time

logger = logging.getLogger("engine")

HAS_UWSGI = False

def int_to_bytes(value):
    return (value).to_bytes(8, byteorder='big')

def int_from_bytes(args, name):
    return int.from_bytes(args[name].encode(), byteorder='big')

COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @timer(20)
    def timer(args):
        from biostar.engine.models import Job
        # Check for queued jobs.
        jobs = Job.objects.filter(state=Job.QUEUED)
        if jobs:
            # Put the jobs in SPOOLED state so that it does not get respooled.
            Job.objects.filter(state=Job.QUEUED).update(state=Job.SPOOLED)
            for job in jobs:
                logger.info(f"Spooling job id={job.id}")
                execute_job.spool(job_id=job.id)

    #@timer(10)
    def spool_demo(args):
        """
        Spooling demo scheduler. Comment out the decorator to launch.
        """
        global COUNTER
        logger.info(f"SPOOL DEMO SUBMIT {COUNTER}, {COUNTER+1}, {COUNTER+2}")
        spool_demo_task.spool(value=COUNTER)
        spool_demo_task.spool(value=COUNTER + 1)
        spool_demo_task.spool(value=COUNTER + 2)
        COUNTER += 3

    @spool(pass_arguments=True)
    def spool_demo_task(value):
        '''
        Spool demonstration task.
        '''
        N = 5
        logger.info(f"-> JOB START {value} ")
        time.sleep(N)
        logger.info(f"<- JOB END {value}")

    @spool(pass_arguments=True)
    def execute_job(job_id):
        """
        Execute job in spooler.
        """
        logger.info(f"Executing spooled job id={job_id}")
        management.call_command('job', id=job_id)



except ModuleNotFoundError as exc:
    pass

def execute(command, workdir="."):
    proc = subprocess.run(command, cwd=workdir, shell=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = force_text(proc.stdout), force_text(proc.stderr)
    return stdout, stderr

