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
        # Find the json and the template.
        json_data = hjson.loads(job.json_data)
        template = job.template

        # This is the work directory.
        workdir = job.path

        # Override template.
        if use_template:
            template = open(use_template).read()

        # Override json.
        if use_json:
            json_data = hjson.loads(open(use_json).read())

        # Print the json.
        if show_json:
            print(hjson.dumps(json_data, indent=4))
            return

        # Print the template.
        if show_template:
            print(template)
            return

        # Extract the execute commands from the spec.
        execute = json_data.get('execute', {})
        filename = execute.get("filename", "Makefile")
        command = execute.get("command", "make all")

        # Render the script.
        template = Template(template)
        context = Context(json_data)
        script = template.render(context)

        # Show the script.
        if show_script:
            print(f'{script}')
            return

        # Logging should start after the early returns.
        logger.info(f'job id={job.id} started.')

        # Make the output directory
        logger.info(f'job id={job.id} workdir: {workdir}')
        if not os.path.isdir(workdir):
            os.mkdir(workdir)

        # Create the script in the output directory.
        with open(os.path.join(workdir, filename), 'wt') as fp:
            fp.write(script)

        # Show the command that is executed.
        logger.info(f'job id={job.id} executing: (cd {workdir} && {command})')

        # Switch the job state to RUNNING.
        job.state = job.RUNNING
        job.save()

        # Run the command.
        proc = subprocess.run(command, cwd=workdir, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Return code indicates an error.
        if (proc.returncode != 0):
            raise Exception(f"executing: {command}")

        # If we made it this far the job has finished.
        job.state = job.FINISHED
        job.save()

    except Exception as exc:
        # Handle all errors here.
        job.state = job.ERROR
        job.save()
        stderr_log.append(f'{exc}')
        logger.info(f'job id={job.id} error {exc}')

    # Collect the output.
    stdout_log.extend(force_text(proc.stdout).splitlines())
    stderr_log.extend(force_text(proc.stderr).splitlines())

    # For now keep logs in one field. TODO: separate into stdin, stdout
    output_log = stdout_log + stderr_log
    job.log = "\n".join(output_log)
    job.save()
    logger.info(f'job id={job.id} finished, status={job.get_state_display()}')

    # Use -v 2 to see the output of the command.
    if verbosity > 1:
        job = Job.objects.get(id=job.id)
        print ("-" * 40)
        print (job.log)
        print("-" * 40)

def error(msg):
    logger.error(msg)
    sys.exit()

class Command(BaseCommand):
    help = 'Job manager.'

    def add_arguments(self, parser):

        parser.add_argument('--next',
                        action='store_true',
                        default=False,
                        help="Runs the oldest queued job")

        parser.add_argument('--id',
                            type=int,
                            default=0,
                            help="Runs job specified by id.")

        parser.add_argument('--show_script',
                            action='store_true',
                            help="Shows the script.")

        parser.add_argument('--show_json',
                            action='store_true',
                            help="Shows the JSON for the job.")

        parser.add_argument('--show_template',
                            action='store_true',
                            help="Shows the template for the job.")

        parser.add_argument('--use_json',
                            help="Override the JSON with this file.")

        parser.add_argument('--use_template',
                            help="Override the TEMPLATE with this file.")

        parser.add_argument('--queued',
                            action='store_true',
                            help="Show most recent 10 queued.")


    def handle(self, *args, **options):

        jobid = options['id']
        next = options['next']
        queued = options['queued']

        if next:
            job = Job.objects.filter(state=Job.QUEUED).order_by('-id').first()
            if not job:
                error(f'there are no queued jobs')
            run(job, options=options)

        if jobid:
            job = Job.objects.filter(id=jobid).first()
            if not job:
                error(f'job for id={jobid} missing')
            run(job, options=options)

        if queued:
            jobs = Job.objects.filter(state=Job.QUEUED).order_by('-id')[:10]
            for job in jobs:
                print(job.id)
            return


