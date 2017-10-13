from django.core.management.base import BaseCommand
from django.template import Template,Context
from engine.models import Job
import subprocess, os, sys,hjson,logging
from django.utils.text import force_text

logger = logging.getLogger('engine')

CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def run(job,show=False):
    ''''
    takes job object, runs the job and return job status
    '''

    logger.info(f'job {job.id} started')
    json_text = job.json_data
    template = job.makefile_template
    outdir = job.path
    error_log = []
    stdout_log =[]
    stderr_log = []

    try:
        # render makefile.
        json_data = hjson.loads(json_text)
        execute = json_data.get('execute',{})
        commands = execute.get("command", "make all").split()
        filename = execute.get("filename", "Makefile")

        template = Template(template)
        context = Context(json_data)

        mtext = template.render(context)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        logger.info(f'job output dir {outdir}')
        with open(os.path.join(outdir, filename), 'wt') as fp:
            fp.write(mtext)

        if show:
            print(f'\n\n{mtext}\n\n')
            return

        # Run the command.
        job.state = job.RUNNING
        logger.info(f'job commands {commands}')
        process = subprocess.run(commands, cwd=outdir, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, check=True)
        job.state = job.FINISHED

    except subprocess.CalledProcessError as err:
        error_log.append(err.stderr.decode('utf-8'))
        job.state = job.ERROR
        logger.info(f'job error: {err}')

    except Exception as exc:
        error_log.append(str(exc))
        job.state = job.ERROR
        logger.info(f'job exception: {exc}')


    print ( force_text(process.stdout) )
    print(force_text(process.stderr))
    #stdout_log.append(process.stdout.read())
    #stderr_log.append(process.stderr.read())
    job.log = "\n".join(stdout_log + error_log)
    job.save()
    logger.info(f'job {job.id} completed')

    return job


class Command(BaseCommand):
    help = 'Run jobs that are queued.'

    def add_arguments(self, parser):

        # positional arguments.
        #parser.add_argument('limit',type=int,default=1,help="Enter the number of jobs to run. Default is 1.")

        # Named (optional) arguments
        parser.add_argument('--run',
                            type =int,
                            default=0,
                            help="Runs job specified by id.")

        parser.add_argument('--show',
                            action='store_true',
                            help="")
        '''
        parser.add_argument('--limit',
                            dest='limit',
                            type =int,
                            default=1, help="Enter the number of jobs to run.")

        parser.add_argument('--jobid',
                             dest='jobid',
                             help="Specifies job id.")

        parser.add_argument('--queued',
                             dest='queued',
                             action='store_true',
                             help="List ten most recent queued jobs.")
        parser.add_argument('--template',
                            dest='template',
                            action='store_true',
                            help="Show template.")
        parser.add_argument('--spec',
                            dest='spec',
                            action='store_true',
                            help="Show analysis spec.")
        '''
    def handle(self, *args, **options):

        #limit = options['limit']
        #jobid = options['jobid']
        #queued = options['queued']
        #spec = options['spec']
        #template = options['template']
        runid = options['run']
        show = options['show']

        if run:
            job = Job.objects.filter(id=runid).first()
            if not job:
                logger.error(f'job for id={runid} missing')
                sys.exit()

            run(job,show=show)
            return
         #jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:limit]

        '''
        if queued:
            jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:10]
            for job in jobs:
                print(job.id)
            return

        if not (jobid or limit) and run:
            print("command requires --jobid or --limit")
            return

        if not jobid and (template or spec):
            print ("command requires --jobid")
            return

        if jobid:
            job = Job.objects.get(id=jobid)

            if template:
                print(job.makefile_template)
                return

            if spec:
                print(job.json_data)
                return
            if run:
                print("runs job")
                run(job)
                return

        if limit:
            jobs = Job.objects.filter(state=Job.QUEUED).order_by("id")[:limit]
            for job in jobs:
                #run(job)
                print("running job")
            return

        '''
