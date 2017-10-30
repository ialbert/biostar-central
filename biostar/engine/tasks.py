'''
This is active only when deployed via UWSGI
'''

import logging, time
from django.core import management

logger = logging.getLogger("engine")

HAS_UWSGI = False

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @timer(5)
    def execute_timer(args):
        from engine.models import Job
        # logger.info(f"executing timer with {args}")
        # It is faster to check here
        first = Job.objects.filter(state=Job.QUEUED).first()
        if first:
            management.call_command('job', next=True)

    @spool
    def execute_job(args):
        job_id = int.from_bytes(args["job_id"].encode(), byteorder='big')
        management.call_command('job', id=job_id)


except ModuleNotFoundError as exc:
    # logger.warning("uwsgidecorators not found, tasks not enabled!")
    pass

