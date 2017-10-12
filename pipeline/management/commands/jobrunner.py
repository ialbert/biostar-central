from django.core.management.base import BaseCommand
from django.template import Template,Context
from engine.models import Job
from pipeline import render
import subprocess, os, sys,hjson


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
        # render makefile.
        spec = hjson.loads(spec)
        template = Template(template)
        context = Context(spec)
        mtext = template.render(context)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        with open(os.path.join(outdir, "Makefile"), 'wt') as fp:
            fp.write(mtext)

        # run makefile.
        process = subprocess.run(['make', 'all'], cwd=outdir, stderr=subprocess.PIPE, check=True)
        return_code = process.returncode

    except subprocess.CalledProcessError as err:

        job.state = job.ERROR
        errorlog.append(err.stderr.decode('utf-8'))
        if len(errorlog) > 100:
            sys.exit(1)

    finally:
        job.log = "\n".join(errorlog)
        with open(os.path.join(outdir, "run_log.txt"), 'wt') as fp:
            fp.write(job.log)
        if job.state != job.ERROR:
            job.state = job.FINISHED
        print(job.state)
        job.save()
    return job


class Command(BaseCommand):
    help = 'Run jobs that are queued.'

    def add_arguments(self, parser):

        # positional arguments.
        parser.add_argument('limit',type=int,default=1,help="Enter the number of jobs to run. Default is 1.")

        # Named (optional) arguments
        parser.add_argument('--show_queued',
                            action='store_true',
                             dest='show_queued',
                            default=False, help="List tem most recent queued jobs.")
        parser.add_argument('--show_make',
                            action='store_true',
                             dest='show_make',
                            default=False, help="Show makefile of the queued jobs.")
        parser.add_argument('--show_spec',
                            action='store_true',
                             dest='show_spec',
                            default=False, help="Show analysis spec of the queued jobs.")

    def handle(self, *args, **options):

        limit = options['limit']
        jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:limit]

        if options['show_make']:
            for job in jobs:
                print("Makefile for job {0}".format(job.id))
                print(job.makefile_template)
            sys.exit(1)

        if options['show_spec']:
            for job in jobs:
                print("Makefile for job {0}".format(job.id))
                print(job.json_data)
            sys.exit(1)

        if options['show_queued']:
            jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:10]
            for job in jobs:
                print(job.id)
            sys.exit(1)

        for job in jobs:
            if job.id ==1:
                run(job)






