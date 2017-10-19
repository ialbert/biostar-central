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

    @timer(60)
    def execute_timer():

        # Run oldest queued job every minute
        from django.core import management
        from .models import Job

        job = Job.objects.filter(state=Job.QUEUED).first()

        management.call_command('job', id=job.id)

    @spool
    def execute_spooler(args):

        import django
        # duh lol
        django.setup()

        jobid = int.from_bytes(args["job_id"].encode(), byteorder='big')

        print(jobid)


except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
