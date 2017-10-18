'''
This is active only when deployed via UWSGI
'''
from django.conf import settings



import logging, time

logger = logging.getLogger('engine')

HAS_UWSGI = False

try:
    from uwsgidecorators import spool

    HAS_UWSGI = True

    @spool
    def execute(args):
        from django.core import management
        #return

        job_id = int.from_bytes(args["jobid"].encode(), byteorder="big")

        #job = args['job']
        #print(management.get_commands())
        print(job_id)
        #management.call_command('job', id=job_id)

except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
