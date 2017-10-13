from django.core.management.base import BaseCommand
from django.template import Template,Context
from engine.models import Job
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
        job.state = job.RUNNING
        job.save()
        process = subprocess.run(['make', 'all'], cwd=outdir, stderr=subprocess.PIPE, check=True)
        return_code = process.returncode

    except subprocess.CalledProcessError as err:

        job.state = job.ERROR
        errorlog.append(err.stderr.decode('utf-8'))
        if len(errorlog) > 100:
            sys.exit(1)

    finally:
        if errorlog:
            job.log = "\n".join(errorlog)
        else:
            job.log = "Analysis completed sucessfully. Results are in results folder. "
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
        #parser.add_argument('limit',type=int,default=1,help="Enter the number of jobs to run. Default is 1.")

        # Named (optional) arguments
        parser.add_argument('--run',
                            action='store_true',
                             dest='run',
                             default =False,
                             help="Runs job. Job should be specified by jobid or limit.")

        parser.add_argument('--limit',
                            action='store_true',
                             dest='limit',
                            default=1, help="Enter the number of jobs to run.")
        parser.add_argument('--jobid',
                            action='store_true',
                             dest='jobid',
                             help="Specifies job id.")

        parser.add_argument('--queued',
                            action='store_true',
                             dest='queued',
                             help="List ten most recent queued jobs.")
        parser.add_argument('--template',
                            action='store_true',
                             dest='template',
                            help="Show template.")
        parser.add_argument('--spec',
                            action='store_true',
                             dest='pec',
                            help="Show analysis spec.")

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
                print("Specs for job {0}".format(job.id))
                print(job.json_data)
            sys.exit(1)

        if options['show_queued']:
            jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:10]
            for job in jobs:
                print(job.id)
            sys.exit(1)

        for job in jobs:
            run(job)






