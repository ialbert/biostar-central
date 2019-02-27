import logging
import hjson
import os
import io
from urllib.parse import urljoin
import requests
import subprocess
import sys
from functools import partial

from django.utils.encoding import force_text
from django.template import Template, Context
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import Q
from django.shortcuts import reverse
from django.utils import timezone

from biostar.emailer.auth import notify
from biostar.engine.models import Analysis, Project, Data, Job
from biostar.engine.api import change_image, get_thumbnail
from biostar.engine import auth
from biostar.accounts.models import User

logger = logging.getLogger('engine')

# Override the logger.
logger.setLevel(logging.INFO)


class Bunch():
    def __init__(self, **kwargs):
        self.value = ''
        self.name = self.summary = ''
        self.help = self.type = self.link = ''
        self.__dict__.update(kwargs)


def build_api_url(root_url, uid=None, view="recipe_api_list", api_key=None):

    url = reverse(view, kwargs=dict(uid=uid)) if uid else reverse(view)
    #TODO Put together params in diffrent way
    full_url = urljoin(root_url, url) + f"?k={api_key}"
    return full_url


def get_json_text(source, target_file=""):

    if os.path.exists(target_file):
        target = hjson.loads(open(target_file, "r").read())
    else:
        target = {}
    # Copy source data into target without overwriting target
    for key in source:
        target[key] = source[key]

    return hjson.dumps(target)


def data_from_json(root, json_data, pid):
    project = Project.objects.get_all(uid=pid).first()

    # The data field is empty.
    if not json_data:
        logger.error(f"JSON file does not have a valid data field")
        return

    # The datalist is built from the json.
    data_list = [Bunch(**row) for row in json_data]

    # Add each collected datatype.
    for bunch in reversed(data_list):
        # This the path to the data.
        path = bunch.value

        # Makes the path relative if necessary.
        path = path if path.startswith("/") else os.path.join(root, path)

        # Create the data.
        auth.create_data(project=project, path=path, type=bunch.type,
                         name=bunch.name, text=bunch.help)


def get_recipes(pid, root_url, api_key):

    if root_url:
        json_url = build_api_url(root_url=root_url, api_key=api_key, view="project_api_info", uid=pid)
        json_text = requests.get(url=json_url, params=dict(k=api_key)).content
        json_data = hjson.loads(json_text)
        # Get the recipes
        recipe_uids = json_data.get("recipes", [])
    else:
        recipes = Analysis.objects.filter(project__uid=pid)
        recipe_uids = recipes.values_list("uid", flat=True)

    if not recipe_uids:
        logger.error(f"No recipes found for pid={pid}")

    return recipe_uids


def generate_fnames(json):

    setting_dict = json.get("settings", {})
    name = setting_dict.get('name')

    name = '_'.join(name.split())
    uid = setting_dict.get("id")
    base = f"{name}-{uid}"

    image = base + '.png'
    template = base + '.sh'
    hjson = base + '.hjson'

    return image, hjson, template


def open_file(absfname, mode="r"):

    if not os.path.exists(absfname):
        logger.error(f"{absfname} does not exist.")
        sys.exit()

    return open(absfname, mode)


def get_project(pid, create=False, privacy=Project.PRIVATE):
    project = Project.objects.get_all(uid=pid).first()
    if not project:
        if create:
            user = User.objects.filter(is_staff=True).first()
            project = auth.create_project(user=user, name="Project Name", uid=pid, privacy=privacy)
        else:
            logger.error(f"*** Project id {pid}: does not exist.")
            sys.exit()
    return project


def get_recipe(rid, pid=None, create=False):
    recipe = Analysis.objects.get_all(uid=rid).first()
    if not recipe:
        if create:
            project = get_project(pid=pid)
            recipe = auth.create_analysis(project=project, json_text="", template="", uid=rid, name="Recipe Name")
        else:
            logger.error(f"*** Recipe id {rid}: does not exist.")
            sys.exit()

    return recipe


def get_response(root_url, view, uid, api_key=""):
    # Send GET request to view and return a response.
    response = requests.get(url=build_api_url(root_url=root_url, api_key=api_key, uid=uid, view=view))
    if response.status_code != 200:
        logger.error(f"*** Error with remote host with id={uid}: {response.text}")
        sys.exit()
    return response


def put_response(root_url, view, uid, stream, api_key=""):
    # Send PUT request to view with data in stream send as a file.
    full_url = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
    response = requests.put(url=full_url, files=dict(file=stream), data=dict(k=api_key))
    if response.status_code != 200:
        logger.error(f"*** Error on remote host: {response.text}")
        sys.exit()
    return response


def start_jobs(recipe, pid):

    project = get_project(pid=pid)
    missing_name = ''
    for key, obj in recipe.json_data.items():
        # When creating a job automatically for data in projects
        # it will try to match the value of the parameter to the data name.
        if obj.get("source") != "PROJECT":
            continue
        name = obj.get('value', '')
        data = Data.objects.filter(project=project, name=name).first()
        if not data:
            missing_name = name
            break
        data.fill_dict(obj)
    if missing_name:
        logger.error(f"Job not created! Missing data:{missing_name} in analysis:{recipe.name}")
    else:
        auth.create_job(analysis=recipe, json_data=recipe.json_data)

    return


def push_recipe(root_dir,  json_file, api_key="", root_url=None, url_from_json=False, jobs=False):
    """
        Push recipe into api/database from a json file.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """

    # Get conf from json file
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    source = abspath(json_file)
    json_data = hjson.loads(open(source, "r").read())
    image_name, template_name, _ = generate_fnames(json=json_data)

    # Get appropriate parameters from conf
    conf = json_data.get("settings", {})
    rid = conf.get("recipe_uid")
    proj_uid = conf.get("project_uid")
    image = conf.get("image") or image_name
    template = conf.get("template") or template_name
    url = conf.get("url") if url_from_json else root_url
    image_stream = open_file(abspath(image), "rb")
    template_stream = open_file(abspath(template), "r")

    if url:
        # Send PUT request to image, json, and template API urls.
        push = lambda view, stream: put_response(root_url=url, view=view, uid=rid, stream=stream, api_key=api_key)
        push(stream=image_stream, view="recipe_api_image")
        push(stream=template_stream, view="recipe_api_template")
        push(stream=open(source, "r"), view="recipe_api_json")
    else:
        recipe = get_recipe(rid=rid, create=True, pid=proj_uid)
        recipe.template = template_stream.read()
        recipe.json_text = hjson.dumps(json_data)
        recipe.name = json_data.get("settings", {}).get("name", recipe.name)
        recipe.text = json_data.get("settings", {}).get("help", recipe.text)
        recipe.image.save(name=image, content=image_stream)
        recipe.save()
        if jobs:
            start_jobs(recipe=recipe, pid=proj_uid)

    print(f"*** Pushed recipe id={rid} from:{source}")

    return


def push_project(root_dir, json_file, root_url=None, api_key="", add_data=False, data_root="", url_from_json=False):
    """
    Load projects from root_dir into remote host or local database
    """
    pmap = {"private": Project.PRIVATE, "public": Project.PUBLIC}
    source = os.path.abspath(os.path.join(root_dir, json_file))
    json_data = hjson.loads(open(source, "r").read())
    image_name, _, _ = generate_fnames(json=json_data)

    conf = json_data.get("settings", {})
    image = conf.get('image') or image_name
    image = os.path.abspath(os.path.join(root_dir, image))
    url = conf.get("url") if url_from_json else root_url
    uid = conf.get("uid")
    privacy = conf.get("privacy", "").lower() or "private"
    privacy = pmap.get(privacy, Project.PRIVATE)
    json_stream = open_file(source)
    image_stream = open_file(image, "rb")

    if url:
        # Send PUT request to image, json, and template API urls.
        push = lambda view, stream: put_response(root_url=url, view=view, uid=uid, stream=stream, api_key=api_key)
        push(stream=json_stream, view="project_api_info")
        push(stream=image_stream, view="recipe_api_json")
    else:
        project = get_project(pid=uid, create=True, privacy=privacy)
        project.name = conf.get("settings", {}).get("name", project.name)
        project.text = conf.get("settings", {}).get("help", project.text)
        project.image.save(name=image, content=image_stream)
        project.save()
        if add_data:
            data = json_data.get("data", [])
            data_from_json(root=data_root, pid=uid, json_data=data)

    print(f"*** Loaded project id=({uid}) from:{source}.")


def pull_recipe(root_dir, rid, root_url=None, api_key="", save=False):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes uid list from API or database.
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    get = lambda view: get_response(root_url=root_url, uid=rid, api_key=api_key, view=view)
    if root_url:
        image = get(view="recipe_api_image").content
        json_text = get(view="recipe_api_json").content.decode()
        template = get(view="recipe_api_template").content.decode()
        json_data = hjson.loads(json_text)
    else:
        recipe = get_recipe(rid=rid)
        json_data = recipe.json_data
        template = recipe.template
        image = open(recipe.image.path if recipe.image else get_thumbnail(), "rb").read()

    print(f"*** Dumping recipe id: {rid}")
    if save:
        # Make output directory.
        os.makedirs(root_dir, exist_ok=True)
        # Write image, template, and json
        img_fname, json_fname, template_fname = generate_fnames(json_data)
        img_fname, json_fname, template_fname = abspath(img_fname), abspath(json_fname), abspath(template_fname)
        open(img_fname, "wb").write(image)
        open(template_fname, "w").write(template)
        open(json_fname, "w").write(hjson.dumps(json_data))
        print(f"{json_fname}\n{img_fname}\n{template_fname}")
    else:
        return json_data, template, image


def pull_project(pid, root_dir, root_url=None, api_key="", save=False):
    """
    Dump project from remote host or local database into root_dir
    """
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    if root_url:
        # Get data from remote url.
        json = get_response(root_url=root_url, api_key=api_key, uid=pid, view="project_api_info").content.decode()
        image = get_response(root_url=root_url, api_key=api_key, uid=pid, view="project_api_image").content
        json = hjson.loads(json)
    else:
        # Get project from database
        project = get_project(pid=pid)
        json = project.json_data
        image = open(project.image.path if project.image else get_thumbnail(), "rb").read()

    print(f"*** Dumped project {pid}: {root_dir}.")
    if save:
        os.makedirs(root_dir, exist_ok=True)
        img_fname, json_fname, _ = generate_fnames(json)
        img_fname, json_fname = abspath(img_fname), abspath(json_fname)
        # Write image to file
        open(img_fname, "wb").write(image)
        # Write hjson to file.
        open(json_fname, "w").write(hjson.dumps(json))
        print(f"{json_fname}\n{abspath(img_fname)}")
    else:
        return json, image

    return True


def data_loader(path, pid=None, uid=None, update_toc=False, name="Data Name", type="", text=""):
    """
    Load data found in path to database.
    """

    data = Data.objects.get_all(uid=uid).first()
    # Work with existing data.
    if data:
        if update_toc:
            data.make_toc()
            print(f"*** Data id : {uid} table of contents updated.")
        return
    project = Project.objects.get_all(uid=pid).first()
    if not project:
        logger.error(f"Project id: {pid} does not exist.")
        return
    if not path or not os.path.exists(path):
        logger.error(f"--path ({path}) does not exist.")
        return
    # Slightly different course of action on file and directories.
    isdir = os.path.isdir(path)

    # Generate alternate names based on input directory type.
    print(f"*** Project: {project.name} ({project.uid})")
    if isdir:
        print(f"*** Linking directory: {path}")
        altname = os.path.split(path)[0].split(os.sep)[-1]
    else:
        print(f"*** Linking file: {path}")
        altname = os.path.split(path)[-1]

    # Get the text from file
    text = open(text, "r").read() if os.path.exists(text) else ""

    # Select the name.
    name = name or altname

    # Create the data.
    data = auth.create_data(project=project, path=path, type=type, name=name, text=text)
    print(f"*** Created data: name={name}, id={data.uid}")

    return data


def get_json_files(root_dir, json_fname=None):
    """
    Return all existing .hjson or .json files in a directory
    """

    is_json = lambda fname: fname.endswith(".hjson") or fname.endswith(".json")
    json_files = [fname.name for fname in os.scandir(root_dir) if fname.is_file() and is_json(fname.name)]
    if json_fname:
        # Return one json file if provided.
        json_files = [json_fname]

    # Filter for files that exist.
    abspath = lambda p: os.path.abspath(os.path.join(root_dir, p))
    recipe_jsons = list(filter(lambda p: os.path.exists(abspath(p)) , json_files))
    return recipe_jsons


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
            print(hjson.dumps(json_data, indent=4))
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
        script_name = execute.get("filename", "recipe.sh")

        # Make the log directory that stores sdout, stderr.
        LOG_DIR = 'runlog'
        log_dir = os.path.join(work_dir, f"{LOG_DIR}")
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # Runtime information will be saved in the log files.
        json_fname = f"{log_dir}/input.json"
        stdout_fname = f"{log_dir}/stdout.txt"
        stderr_fname = f"{log_dir}/stderr.txt"

        # Build the command line
        command = execute.get("command", "bash recipe.sh")

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
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        # Create the script in the output directory.
        with open(os.path.join(work_dir, script_name), 'wt') as fp:
            fp.write(script)

        # Create a file that stores the json data for reference.
        with open(json_fname, 'wt') as fp:
            fp.write(hjson.dumps(json_data, indent=4))

        # Initial create each of the stdout, stderr file placeholders.
        for path in [stdout_fname, stderr_fname]:
            with open(path, 'wt') as fp:
                pass

        # Show the command that is executed.
        logger.info(f'Job id={job.id} executing: {full_command}')

        # Job must be authorized to run.
        if job.security != Job.AUTHORIZED:
            raise Exception(f"Job security error: {job.get_security_display()}")

        # Switch the job state to RUNNING and save the script field.
        Job.objects.filter(pk=job.pk).update(state=Job.RUNNING,
                                             start_date=timezone.now(),
                                             script=script)
        # Run the command.
        proc = subprocess.run(command, cwd=work_dir, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Raise an error if returncode is anything but 0.
        proc.check_returncode()

        # If we made it this far the job has finished.
        logger.info(f"uid={job.uid}, name={job.name}")
        Job.objects.filter(pk=job.pk).update(state=Job.COMPLETED)

    except Exception as exc:
        # Handle all errors here.
        Job.objects.filter(pk=job.pk).update(state=Job.ERROR)
        stderr_log.append(f'{exc}')
        logger.error(f'job id={job.pk} error {exc}')

    # Collect the output.
    if proc:
        stdout_log.extend(force_text(proc.stdout).splitlines())
        stderr_log.extend(force_text(proc.stderr).splitlines())

    # Save the logs and end time
    Job.objects.filter(pk=job.pk).update(end_date=timezone.now(),
                                         stdout_log="\n".join(stdout_log),
                                         stderr_log="\n".join(stderr_log))

    # Reselect the job to get refresh fields.
    job = Job.objects.filter(pk=job.pk).first()

    # Create a log script in the output directory as well.
    with open(stdout_fname, 'wt') as fp:
        fp.write(job.stdout_log)

    # Create a log script in the output directory as well.
    with open(stderr_fname, 'wt') as fp:
        fp.write(job.stderr_log)

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
        notify(template_name="emailer/job_finished.html", email_list=[job.owner.email], send=True,
               extra_context=context)


def list_obj(mtype="project"):

    mtype_map = dict(project=Project, recipe=Analysis, data=Data, job=Job)

    objs = mtype_map.get(mtype).objects.get_all().order_by("id")[:100]
    for instance in objs:
        print(f'{instance.id}\t{instance.uid}\t{instance.project.uid}\t{instance.name}')

    return


class Command(BaseCommand):
    help = 'Dump and load items using api.'

    def load_data(self, options):
        root_dir = options.get("dir")
        pid = options.get("pid")
        did = options.get("did")
        path = options.get("path") or ""
        name = options.get("name")
        type = options.get("type")
        update_toc = options.get("update_toc")

        if not (pid or did):
            self.stdout.write(self.style.NOTICE("[error] --pid or --did need to be set."))
            self.run_from_argv(sys.argv + ["--help"])
            sys.exit()
        path = path or os.path.abspath(root_dir)
        data = data_loader(pid=pid, path=path, uid=did, update_toc=update_toc, name=name, type=type)
        msg = f"{data.name} loaded into database."
        self.stdout.write(msg=self.style.SUCCESS(msg))
        return

    def manage_job(self, **options):
        jobid = options.get('id')
        jobuid = options.get('uid')
        next = options.get('next')
        queued = options.get('list')

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
                logger.info(f'job for id={jobid}/uid={jobuid} missing')
            else:
                run(job.first(), options=options)
            return

        if queued:
            jobs = Job.objects.get_all().order_by('id')[:100]
            for job in jobs:
                print(f'{job.id}\t{job.get_state_display()}\t{job.name}')
            return
        return

    def manage_push(self, **options):
        subcommand = sys.argv[2] if len(sys.argv) > 2 else None
        push = subcommand == "push"
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir")
        rid = options.get("rid")
        pid = options.get("pid") or ""

        data = options.get("data")
        did = options.get("did")

        create_job = options.get("jobs")
        json_file = options.get("json")
        data_root = options.get("data_root")
        add_data = options.get("data_from_json")
        url_from_json = options.get("url_from_json")

        # Handle loading data
        if data or did:
            self.load_data(options=options)
            return

        # Get root dir from dir name of json file.
        if json_file:
            full_path = os.path.abspath(json_file)
            root_dir = root_dir or os.path.dirname(full_path)
            json_file = os.path.basename(full_path)

        # Require api key when pushing to remote url
        if ((root_url or url_from_json) and push) and not api_key:
            sys.argv.append("--help")
            self.stdout.write(self.style.NOTICE("[error] --key is required when loading data to remote site."))
            self.run_from_argv(sys.argv)
            sys.exit()

        root_dir = root_dir or os.getcwd()
        if not os.path.exists(root_dir):
            logger.error(f"*** Directory: {root_dir} does not exist.")
            sys.exit()

        # Get json files from the root dir.
        json_files = get_json_files(root_dir=root_dir, json_fname=json_file)
        print(f"Pushing files from --dir {root_dir}")
        for fname in json_files:
            json_text = open(os.path.abspath(os.path.join(root_dir, fname)), "r").read()
            recipe_uid = hjson.loads(json_text).get("settings", {}).get("recipe_uid")
            project_uid = hjson.loads(json_text).get("settings", {}).get("project_uid")

            # Skip pushing when rid/pid in json != --rid/--pid given
            if (rid and recipe_uid != rid) or (pid and project_uid != pid):
                continue
            if recipe_uid:
                push_recipe(root_dir=root_dir, root_url=root_url, api_key=api_key, json_file=fname,
                            jobs=create_job, url_from_json=url_from_json)
            else:
                push_project(root_dir=root_dir, root_url=root_url, api_key=api_key, add_data=add_data,
                             data_root=data_root, url_from_json=url_from_json, json_file=fname)

        logger.info(f"{len(json_files)} pushed into {'url' if (url_from_json or root_url) else 'database'}.")

        return

    def manage_pull(self, **options):

        pull_recipes = options.get("recipes")
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir") or os.getcwd()
        rid = options.get("rid")
        pid = options.get("pid")

        if (not (pid or rid)) or len(sys.argv) == 3:
            sys.argv.append("--help")
            self.stdout.write(self.style.NOTICE("--pid or --rid is required."))
            self.run_from_argv(sys.argv)
            sys.exit()

        print(f"Writing into directory: {root_dir}.")
        if rid:
            pull_recipe(root_dir=root_dir, root_url=root_url, api_key=api_key, rid=rid, save=True)
            logger.info(f"Recipe id {rid} dumped into {root_dir}.")
            return

        if pull_recipes:
            # Get recipes belonging to project --pid
            recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key)
            # Get multiple recipes belonging to project --pid
            for uid in recipes:
                pull_recipe(root_dir=root_dir, root_url=root_url, api_key=api_key, rid=uid, save=True)
                logger.info(f"Recipe id {uid} dumped into {root_dir}.")
            return

        pull_project(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key, save=True)
        logger.info(f"Project id: {pid} dumped into {root_dir}.")
        return

    def add_push_commands(self, parser):

        parser.add_argument('-u', "--url_from_json", action="store_true",
                            help="""Extract url from conf file instead of --url.""")
        parser.add_argument('-d', "--data", action="store_true",
                            help="""Load --dir or --path as a data object of --pid to local database.""")
        self.add_api_commands(parser=parser)

        parser.add_argument("--data_from_json", action='store_true', help="Add data found in --json to --pid.")
        parser.add_argument("--jobs", action='store_true', help="Also creates a queued job for the recipe")
        parser.add_argument('--rid', type=str, default="", help="Recipe uid to load.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to load from or dump to.")
        parser.add_argument("--did", type=str, help="Data uid to load or update.")

        parser.add_argument('--dir', default='', help="Base directory to store/load from.")
        parser.add_argument("--path", type=str, help="Path to data.")
        parser.add_argument('--text', default='', help="A file containing the description of the data")
        parser.add_argument('--name', default='', help="Sets the name of the data")
        parser.add_argument('--type', default='data', help="Sets the type of the data")
        parser.add_argument("--update_toc", action="store_true", help="Update table of contents for data --uid.")

        parser.add_argument("--data_root", default="",
                            help="Root directory to data found in conf file when loading project.")
        parser.add_argument('--json', default='', help="""JSON file path relative to --dir to get conf from.""")
        return

    def add_pull_commands(self, parser):
        parser.add_argument('-r', "--recipes", action="store_true",
                            help="""Pull recipes of --pid""")
        self.add_api_commands(parser=parser)

        parser.add_argument('--rid', type=str, default="", help="Recipe uid to dump.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to dump.")
        parser.add_argument('--dir', default='', help="Directory to store in.")
        return

    def add_api_commands(self, parser):
        """Add default api commands to parser"""
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")

        return

    def add_job_commands(self, parser):

        parser.add_argument('--next', action='store_true', default=False, help="Runs the oldest queued job")
        parser.add_argument('--id', type=int, default=0, help="Runs job specified by id.")
        parser.add_argument('--jid', type=str, default='', help="Runs job specified by uid.")
        parser.add_argument('--show_script', action='store_true', help="Shows the script.")
        parser.add_argument('--show_json', action='store_true', help="Shows the JSON for the job.")
        parser.add_argument('--show_template', action='store_true', help="Shows the template for the job.")
        parser.add_argument('--show_command', action='store_true', help="Shows the command executed for the job.")
        parser.add_argument('--use_json', help="Override the JSON with this file.")
        parser.add_argument('--use_template', help="Override the TEMPLATE with this file.")
        parser.add_argument('--list', action='store_true', help="Show a job list")

    def add_arguments(self, parser):

        subparsers = parser.add_subparsers()

        load_parser = subparsers.add_parser("push", help="""
                                                    Project, Data and Recipe push manager.
                                                    Project:  Create or update project in remote host or database.
                                                    Data:     Create or update data to project --pid in database. 
                                                    Recipe:   Create or update recipe in remote host or database. 
                                                    ."""
                                            )
        self.add_push_commands(parser=load_parser)

        dump_parser = subparsers.add_parser("pull", help="""
                                                    Project and Recipe Job dumper manager.
                                                    Project  : Dump project from remote host or database
                                                    Recipe:  : Dump Recipe from remote host or database.
                                                    """)
        self.add_pull_commands(parser=dump_parser)

        job_parser = subparsers.add_parser("job", help="Job manager.")
        self.add_job_commands(parser=job_parser)

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None

        listing = options.get("list")

        if len(sys.argv) <= 3:
            sys.argv.append("--help")
            self.run_from_argv(sys.argv)
            sys.exit()

        if listing:
            list_obj(mtype=subcommand)
            return

        if subcommand == "push":
            self.manage_push(**options)
            return

        if subcommand == "pull":
            self.manage_pull(**options)
            return

        if subcommand == "job":
            self.manage_job(**options)
            return
