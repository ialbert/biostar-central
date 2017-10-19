'''
This is active only when deployed via UWSGI
'''

from django.conf import settings
import logging, time

logger = logging.getLogger('engine')

HAS_UWSGI = False

try:
    from uwsgidecorators import *

    HAS_UWSGI = True

    @timer(3)
    def execute_timer(job_id):

        from django.core import management

        management.call_command('job', id=job_id)

    @spool
    def execute_spooler(args):

        jobid = int.from_bytes(args["job_id"].encode(), byteorder='big')

        print(jobid)


except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
