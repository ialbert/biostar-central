import toml as hjson
import time
import os, logging, subprocess, pprint

from django.conf import settings
from django.core.management.base import BaseCommand
from django.template import Template, Context
from django.utils.encoding import force_text

from biostar.recipes.models import Job
from biostar.recipes import auth
from django.conf import settings
from django.utils import timezone
from biostar.emailer.tasks import send_email

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.DEBUG)

CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def finalize_job(job, data):
    """
    Insert file from data dict into database.

    # Only file is required, it can be a directory as well.
    [[settings.create]]
    file = runlog
    uid =
    name
    text =
    type =

     # List of files
    {'settings': {'create': [ {'file':'foo', }, {'file':'bar'} ]  }

    """

    # Get files intended for insertion.
    create = data.get("settings", {}).get("create", [])
    root = job.path
    project = job.project

    for create_using in create:

        # Create full path using the file.
        fname = create_using.get("file", '')
        fname = fname.strip()
        fullpath = os.path.abspath(os.path.join(root, fname))

        # Ensure the file path is valid.
        valid = len(fname) and os.path.exists(fullpath) and fullpath.startswith(root)

        if not valid:
            raise FileNotFoundError(f"File: {fname} does not exist.")

        # Add the full path and project to data dict.
        create_using.update({"file": fullpath, "project": project})

        auth.get_or_create(**create_using)


def create_logs(job):

    work_dir = job.path

    # Runtime information will be saved in the log files.
    stdout_fname = os.path.join(work_dir, settings.JOB_STDOUT)
    stderr_fname = os.path.join(work_dir, settings.JOB_STDERR)

    # Initial create each of the stdout, stderr file placeholders.
    for path in [stdout_fname, stderr_fname]:
        dirname = os.path.dirname(path)
        os.makedirs(dirname, exist_ok=True)
        with open(path, 'wt') as fp:
            pass

    return stdout_fname, stderr_fname


def run(job, options={}):
    """
    Runs a job
    """
    # Options that cause early termination.
    show_json = options.get('show_json')
    show_template = options.get('show_template')
    show_script = options.get('show_script')
    show_command = options.get('show_command')
    use_template = options.get('use_template')
    use_json = options.get('use_json')
    verbosity = options.get('verbosity', 0)

    # Create log directories and files
    stdout_fname, stderr_fname = create_logs(job)

    try:
        # Find the json and the template.
        json_data = hjson.loads(job.json_text)
        template = job.template

        # This is the work directory.
        work_dir = job.path

        # The bade URL of the site.
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
            print(hjson.dumps(json_data))
            return

        # Print the template.
        if show_template:
            print(template)
            return

        # Extract the execute commands from the spec.
        settings_dict = json_data.get("settings", {})

        # Specifies the command that gets executed.
        execute = settings_dict.get('execute', {})

        # The name of the file that contain the commands.
        script_name = execute.get("script_name", "recipe.sh")

        # Runtime information will be saved in the log files.
        json_fname = os.path.join(job.path, f"{settings.JOB_LOGDIR}", "input.json")

        # Build the command line
        command = execute.get("command", f"bash {script_name}")

        # The commands can be substituted as well.
        context = Context(json_data)
        command_template = Template(command)
        command = command_template.render(context)

        # This is the full command that will be executed.
        full_command = f'(cd {work_dir} && {command})'
        if show_command:
            print(full_command)
            return

        # Script template.
        context = Context(json_data)
        script_template = Template(template)
        script = script_template.render(context)

        # Show the script.
        if show_script:
            print(f'{script}')
            return

        # Logging should start after the early returns.
        logger.info(f'Job id={job.id} name={job.name}')

        # Make the output directory
        logger.info(f'Job id={job.id} work_dir: {work_dir}')
        os.makedirs(work_dir, exist_ok=True)

        # Create the script in the output directory.
        with open(os.path.join(work_dir, script_name), 'wt') as fp:
            fp.write(script)

        # Create a file that stores the json data for reference.
        with open(json_fname, 'wt') as fp:
            fp.write(hjson.dumps(json_data))

        # Show the command that is executed.
        logger.info(f'Job id={job.id} executing: {full_command}')

        # Job must be authorized to run.
        if job.security != Job.AUTHORIZED:
            raise Exception(f"Job security error: {job.get_security_display()}. Recipe security : {job.analysis.get_security_display()}")

        # Switch the job state to RUNNING and save the script field.
        Job.objects.filter(pk=job.pk).update(state=Job.RUNNING,
                                             start_date=timezone.now(),
                                             script=script)
        # Run the command.
        proc = subprocess.run(command, cwd=work_dir, shell=True,
                              stdout=open(stdout_fname, "w"),
                              stderr=open(stderr_fname, "w"))

        # Raise an error if returncode is anything but 0.
        proc.check_returncode()

        # Perform tasks at job finalization
        finalize_job(data=json_data, job=job)

        # If we made it this far the job has finished.
        logger.info(f"uid={job.uid}, name={job.name}")
        Job.objects.filter(pk=job.pk).update(state=Job.COMPLETED)

    except Exception as exc:
        # Write error to log file
        open(stderr_fname, "a").write(f"\n{exc}")
        # Handle all errors here.
        Job.objects.filter(pk=job.pk).update(state=Job.ERROR)
        logger.error(f'job id={job.pk} error {exc}')

    stdout_log = open(stdout_fname, "r").read()
    stderr_log = open(stderr_fname, "r").read()
    # Save the logs and end time
    Job.objects.filter(pk=job.pk).update(end_date=timezone.now(),
                                         stdout_log=stdout_log,
                                         stderr_log=stderr_log)

    # Reselect the job to get refresh fields.
    job = Job.objects.filter(pk=job.pk).first()

    # Log job status.
    logger.info(f'Job id={job.id} finished, status={job.get_state_display()}')

    # Use -v 2 to see the output of the command.
    if verbosity > 1:
        print("-" * 40)
        print(job.stdout_log)
        print("-" * 40)
        print(job.stderr_log)

    if job.owner.profile.notify:

        context = dict(subject=job.project.name, job=job)

        # Send notification emails
        send_email(template_name="emailer/job_finished.html", recipient_list=[job.owner.email],
                   extra_context=context)


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
        parser.add_argument('--uid',
                            type=str,
                            default='',
                            help="Runs job specified by uid.")

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
        jobuid = options['uid']
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

        if jobid or jobuid:

            job = Job.objects.filter(uid=jobuid) or Job.objects.filter(id=jobid)
            if not job:
                logger.info(f'job for id={jobid}/uid="{jobuid}"" missing')
            else:
                run(job.first(), options=options)
            return

        if queued:
            jobs = Job.objects.all().order_by('id')[:100]
            for job in jobs:
                print(f'{job.id}\t{job.get_state_display()}\t{job.name}')
            return
