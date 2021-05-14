'''
This is active only when deployed via UWSGI
'''

import logging, time, shutil, subprocess
from django.core import management
from biostar.utils.decorators import task, timer

import time

logger = logging.getLogger("engine")


@timer(250)
def send_emails(*args ,**kwargs):
    """
    Sends  queued emails
    """
    try:
        from mailer.engine import send_all
        send_all()
        logger.info("send_all()")

    except Exception as exce:
        logger.error(exce)


# @timer(2)
# def test(*args ,**kwargs):
#     from biostar.recipes.models import Job
#     print(Job.objects.all())
#     print('TESTING ' *10)
#     logger.info('TESTING'*10)
#

@timer(30)
def scheduler(*args, **kwargs):
    try:
        from biostar.recipes.models import Job
        # Check for queued jobs.
        jobs = Job.objects.filter(state=Job.QUEUED)
        if jobs:
            # Put the jobs in SPOOLED state so that it do not get respooled.
            for job in jobs:
                execute_job.spool(job_id=job.id)

            jobs.update(state=Job.SPOOLED)
    except Exception as exc:
        logger.error(exc)


#@timer(10)
# def spool_demo(args):
#     """
#     Spooling demo scheduler. Comment out the decorator to launch.
#     """
#     global COUNTER
#     logger.info(f"SPOOL DEMO SUBMIT {COUNTER}, {COUNTER+1}, {COUNTER+2}")
#     spool_demo_task.spool(value=COUNTER)
#     spool_demo_task.spool(value=COUNTER + 1)
#     spool_demo_task.spool(value=COUNTER + 2)
#     COUNTER += 3
#
#
# @spool(pass_arguments=True)
# def spool_demo_task(value):
#     '''
#     Spool demonstration task.
#     '''
#     N = 5
#     logger.info(f"-> JOB START {value} ")
#     time.sleep(N)
#     logger.info(f"<- JOB END {value}")


@task
def execute_job(job_id):
    """
    Execute job in spooler.
    """
    logger.info(f"Executing spooled job id={job_id}")
    from biostar.recipes.models import Job

    # Spend 3 seconds in queued state.
    time.sleep(3)
    Job.objects.filter(id=job_id).update(state=Job.SPOOLED)
    # Spend 3 seconds in spooled state
    time.sleep(3)

    management.call_command('job', id=job_id)

#
# @timer(2)
# def foo(name, value):
#    print (f"TIMER TEST {name}={value}")
# foo.timer("bar", value=100)
#
