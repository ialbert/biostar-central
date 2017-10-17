'''
This is active only when deployed via UWSGI
'''
from django.conf import settings
from django.core.management import call_command


import logging, time

logger = logging.getLogger('engine')

HAS_UWSGI = False

try:
    from uwsgidecorators import spool

    HAS_UWSGI = True

    @spool
    def execute(args):
        #return

        job_id = int.from_bytes(args["jobid"].encode(), byteorder="big")


        #job = args['job']
        call_command('job', id=job_id)

except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
