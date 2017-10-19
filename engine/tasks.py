'''
This is active only when deployed via UWSGI
'''

import logging, time
from django.core import management

logger = logging.getLogger("biostar.tasks")

HAS_UWSGI = False

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @timer(30)
    def execute_timer(job_id):
        logger.info("executing timer")


    @spool
    def execute_job(args):

        job_id = int.from_bytes(args["job_id"].encode(), byteorder='big')
        management.call_command('job', id=job_id)


except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
