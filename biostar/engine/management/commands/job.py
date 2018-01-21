import hjson
import os, logging, subprocess, pprint

from django.conf import settings
from django.core.management.base import BaseCommand
from django.template import Template, Context
from django.utils.text import force_text

from biostar.engine.models import Job
from django.utils import timezone

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.DEBUG)

CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def run(job, options={}):
    '''
    Runs a job
    '''
    # Options that cause early termination.
    show_json = options.get('show_json')
    show_template = options.get('show_template')
    show_script = options.get('show_script')
    show_command = options.get('show_command')
    use_template = options.get('use_template')
    use_json = options.get('use_json')
    verbosity = options.get('verbosity', 0)

    # Defined in case we bail on errors before setting it.
    script = command = proc = None

    stdout_log = []
    stderr_log = []
    try:
        # Find the json and the template.
        json_data = hjson.loads(job.json_text)
        template = job.template

        # This is the work directory.
        work_dir = job.path

        url_base = f'{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}'
        # Populate extra context
        def extra_context(job):
            extras = dict(
                media_root=settings.MEDIA_ROOT,
                media_url=settings.MEDIA_URL,
                work_dir=work_dir, local_root=settings.LOCAL_ROOT,
                user_id=job.owner.id, user_email=job.owner.email,
                job_id=job.id, job_name=job.name,
                job_url=f'{url_base}{settings.MEDIA_URL}{job.get_url()}'.rstrip("/"),
                project_id=job.project.id, project_name=job.project.name,
                analyis_name=job.analysis.name,
                analysis_id=job.analysis.id,
                domain=settings.SITE_DOMAIN, protocol=settings.PROTOCOL,
            )
            return extras

        # Add the runtime context.
        json_data['runtime'] = extra_context(job)

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
        settings_dict = json_data.get("settings", {})
        execute = settings_dict.get('execute', {})

        pprint.pprint(execute)

        script_name = execute.get("filename", "run.sh")
        json_fname = f"{script_name}.json"
        stdout_fname = f"{script_name}.stdout.txt"
        stderr_fname = f"{script_name}.stderr.txt"

        command = execute.get("command", "bash run.sh")

        # This is the full command that will be executed.
        full_command = f'(cd {work_dir} & {command})'
        if show_command:
            print(full_command)
            return

        template = Template(template)
        context = Context(json_data)
        script = template.render(context)

        # Store the script for the job.
        job.script = script

        # Show the script.
        if show_script:
            print(f'{script}')
            return

        # Logging should start after the early returns.
        logger.info(f'Job id={job.id} name={job.name}')

        # Make the output directory
        logger.info(f'Job id={job.id} work_dir: {work_dir}')
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        # Create the script in the output directory.
        with open(os.path.join(work_dir, script_name), 'wt') as fp:
            fp.write(script)

        # Create a file that stores the json data for reference.
        with open(os.path.join(work_dir, json_fname), 'wt') as fp:
            fp.write(hjson.dumps(json_data, indent=4))

        # Show the command that is executed.
        logger.info(f'Job id={job.id} executing: {full_command}')

        # Job must be authorized to run.
        if job.security != Job.AUTHORIZED:
            raise Exception(f"Job security error: {job.get_security_display()}")

        # Switch the job state to RUNNING.
        job.state = Job.RUNNING
        job.start_date = timezone.now()
        job.save()

        # Run the command.
        proc = subprocess.run(command, cwd=work_dir, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Return code indicates an error.
        if (proc.returncode != 0):
            raise Exception(f"executing: {command}")

        # If we made it this far the job has finished.
        job.state = Job.COMPLETED
        job.save()

    except Exception as exc:
        # Handle all errors here.
        job.state = Job.ERROR
        job.save()
        stderr_log.append(f'{exc}')
        logger.info(f'job id={job.id} error {exc}')

    # Collect the output.
    if proc:
        stdout_log.extend(force_text(proc.stdout).splitlines())
        stderr_log.extend(force_text(proc.stderr).splitlines())

    # Save the logs.
    job.stdout_log = "\n".join(stdout_log)
    job.stderr_log = "\n".join(stderr_log)

    # Set the end time
    job.end_date = timezone.now()
    job.save()

    # Create a log script in the output directory as well.
    with open(os.path.join(work_dir, stdout_fname), 'wt') as fp:
        fp.write(job.stdout_log)

    # Create a log script in the output directory as well.
    with open(os.path.join(work_dir, stderr_fname), 'wt') as fp:
        fp.write(job.stderr_log)

    logger.info(f'Job id={job.id} finished, status={job.get_state_display()}')
    # Use -v 2 to see the output of the command.
    if verbosity > 1:
        job = Job.objects.get(id=job.id)
        print("-" * 40)
        print(job.stdout_log)
        print("-" * 40)
        print(job.stderr_log)


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

        parser.add_argument('--show_command',
                            action='store_true',
                            help="Shows the command executed for the job.")

        parser.add_argument('--use_json',
                            help="Override the JSON with this file.")

        parser.add_argument('--use_template',
                            help="Override the TEMPLATE with this file.")

        parser.add_argument('--list',
                            action='store_true',
                            help="Show a job list")

    def handle(self, *args, **options):

        jobid = options['id']
        next = options['next']
        queued = options['list']

        # This code is also run insider tasks.
        if next:
            job = Job.objects.filter(state=Job.QUEUED).order_by('id').first()
            if not job:
                logger.info(f'there are no queued jobs')
            else:
                run(job, options=options)
            return

        if jobid:
            job = Job.objects.filter(id=jobid).first()
            if not job:
                logger.info(f'job for id={jobid} missing')
            else:
                run(job, options=options)
            return

        if queued:
            jobs = Job.objects.all().order_by('id')[:100]
            for job in jobs:
                print(f'{job.id}\t{job.get_state_display()}\t{job.name}')
            return
