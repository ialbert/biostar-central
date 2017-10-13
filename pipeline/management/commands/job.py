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
    error_log = []
    stdout_log =[]
    stderr_log = []

    try:
        # render makefile.
        spec = hjson.loads(spec)
        execute = spec.get('execute',{})
        commands = execute.get("command", "make all").split()
        filename = execute.get("filename", "Makefile")


        template = Template(template)
        context = Context(spec)
        mtext = template.render(context)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        with open(os.path.join(outdir, filename), 'wt') as fp:
            fp.write(mtext)

        # Run the command.
        job.state = job.RUNNING
        process = subprocess.run(commands, cwd=outdir, stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True)
        job.state = job.FINISHED


    except subprocess.CalledProcessError as err:
        error_log.append(err.stderr.decode('utf-8'))
        job.state = job.ERROR

    except Exception as exc:
        error_log.append(str(exc))
        job.state = job.ERROR

    finally:
        stdout_log.append(process.stdout.read())
        stderr_log.append(process.stderr.read())
        job.log = "\n".join(stdout_log + error_log)
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






