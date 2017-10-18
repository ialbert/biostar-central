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

    # run every 3 seconds
    @timer(3)
    def execute(job_id):

        from .management.commands import job
        from .models import Job
        #from django.core import management


        #job_id = int.from_bytes(args["jobid"].encode(), byteorder="big")

        #job = args['job']
        #print(management.get_commands())
        #print(job_id)
        current_job = Job.objects.filter(id=job_id).first()

        if current_job:
            job.run(current_job)
        #management.call_command('job', id=job_id)

except ImportError as exc:


    logger.warning("uwsgidecorators not found, tasks not enabled!")
