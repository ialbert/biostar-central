from django.core.management.base import BaseCommand
from django.template import Template, Context
from engine.models import Job
import subprocess, os, sys, hjson, logging
from django.utils.text import force_text

logger = logging.getLogger('engine')

CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def run(job, options={}):
    ''''
    Runs a json
    '''

    # Options that cause early termination.
    show_json = options.get('show_json')
    show_template = options.get('show_template')
    show_script = options.get('show_script')
    use_template = options.get('use_template')
    use_json = options.get('use_json')
    verbosity = options.get('verbosity')

    stdout_log = []
    stderr_log = []

    try:
        logger.info(f'job id={job.id} started.')

        # Extract the used data.
        json_data = hjson.loads(job.json_data)
        template = job.makefile_template
        workdir = job.path

        if use_template:
            template = open(use_template).read()

        if use_json:
            json_data = hjson.loads(open(use_json).read())

        if show_json:
            print(hjson.dumps(json_data, indent=4))
            return

        if show_template:
            print(template)
            return

        # Extract the execute commands from the spec.
        execute = json_data.get('execute', {})
        filename = execute.get("filename", "Makefile")
        command = execute.get("command", "make all")

        template = Template(template)
        context = Context(json_data)
        script = template.render(context)

        # Show the script.
        if show_script:
            print(f'{script}')
            return

        # Make the output directory
        logger.info(f'job id={job.id} workdir: {workdir}')
        if not os.path.isdir(workdir):
            os.mkdir(workdir)

        # Create the script in the output directory.
        with open(os.path.join(workdir, filename), 'wt') as fp:
            fp.write(script)

        # Show the command that is executed.
        logger.info(f'job id={job.id} executing: (cd {workdir} && {command})')

        # Update the job state.
        job.state = job.RUNNING
        job.save()

        # Run the command.
        proc = subprocess.run(command, cwd=workdir, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        job.state = job.FINISHED
        job.save()

        # Return code indicates an error.
        if (proc.returncode != 0):
            raise Exception(f"error when executing: {command}")

    except Exception as exc:
        job.state = job.ERROR
        job.save()
        stderr_log.append(f'{exc}')
        logger.info(f'job id={job.id} {exc}')

    # Collect the output.
    stdout_log.extend(force_text(proc.stdout).splitlines())
    stderr_log.extend(force_text(proc.stderr).splitlines())

    # For now keep it in one location.
    output_log = stdout_log + stderr_log
    job.log = "\n".join(output_log)
    job.save()
    logger.info(f'job id={job.id} finished, status={job.get_state_display()}')

    if verbosity > 1:
        job = Job.objects.get(id=job.id)
        print ("-" * 40)
        print (job.log)
        print("-" * 40)





class Command(BaseCommand):
    help = 'Run jobs that are queued.'

    def add_arguments(self, parser):

        # positional arguments.
        # parser.add_argument('limit',type=int,default=1,help="Enter the number of jobs to run. Default is 1.")

        # Named (optional) arguments
        parser.add_argument('--run',
                            type=int,
                            default=0,
                            help="Runs job specified by id.")

        parser.add_argument('--show_script',
                            action='store_true',
                            help="")

        parser.add_argument('--show_json',
                            action='store_true',
                            help="")

        parser.add_argument('--show_template',
                            action='store_true',
                            help="")

        parser.add_argument('--use_json',
                            help="Use the specified file as JSON data.")

        parser.add_argument('--use_template',
                            help="Use the specified TEMPLATE to render the data.")

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

        # limit = options['limit']
        # jobid = options['jobid']
        # queued = options['queued']
        # spec = options['spec']
        # template = options['template']
        runid = options['run']

        if runid:
            job = Job.objects.filter(id=runid).first()

            if not job:
                logger.error(f'job for id={runid} missing')
            run(job, options=options)

            return

        # jobs = Job.objects.filter(state=Job.QUEUED).order_by("-id")[:limit]

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
