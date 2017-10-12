from django.core.management.base import BaseCommand
from django.template import Template,Context
from engine.models import Job
from pipeline import render
import subprocess, os, sys


CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def run(job):
    ''''
    takes job object, runs the job and return job status
    '''

    spec = job.json_data
    template = job.makefile_template
    outdir = job.path
    errorlog = []

    try:

        mtext = render.render_data(spec,template)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        with open(os.path.join(outdir, "Makefile"), 'wt') as fp:
            fp.write(mtext)

        process = subprocess.run(['make', 'all'], cwd=outdir, stderr=subprocess.PIPE, check=True)
        job.status = process.returncode

    except subprocess.CalledProcessError as err:
        print("ERROR!!")
        errorlog.append(err.stderr.decode('utf-8'))
        if len(errorlog) > 100:
            sys.exit(1)
        job.status = err.returncode

    finally:
        print(CURR_DIR)
        job.log = "\n".join(errorlog)
        print("printing errlog")
        print(job.log)
        #job.save()
    return job


class Command(BaseCommand):
    help = 'Run jobs that are queued.'

    def add_arguments(self, parser):
        parser.add_argument('limit', default=1, type=int)

    def handle(self, *args, **options):
        limit = options['limit']
        jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:limit]
        for job in jobs:
            run(job)





