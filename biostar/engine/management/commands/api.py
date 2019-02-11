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


def get_json_data(json_text, override_json=False, outfile=""):

    # Completely override json in outfile with json_text
    if override_json:
        data = hjson.dumps(hjson.loads(json_text))
        return data

    # Only replace name and help in outfile with items in json_text
    if os.path.exists(outfile):
        current_json = hjson.loads(open(outfile, "r").read())
    else:
        current_json = {}

    name = hjson.loads(json_text).get("settings", {}).get("name", "Name")
    text = hjson.loads(json_text).get("settings", {}).get("help", "Text")

    if current_json.get("settings"):
        current_json["settings"]["name"] = name
        current_json["settings"]["help"] = text
    else:
        current_json["settings"] = {'name': name, 'help': text}

    return hjson.dumps(current_json)


def remote_upload(stream, root_url, uid, api_key, view):
    """
    Upload data found in stream to root_url.
    Currently uses PUT requests
    """

    payload = dict(k=api_key)
    # Build api url then send PUT request.
    full_url = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
    response = requests.put(url=full_url, files=dict(file=stream), data=payload)
    if response.status_code == 404:
        logger.error(f"*** Object id : {uid} does not exist on remote host.")
        sys.exit()

    return response


def remote_download(root_url, api_key, view, uid, is_image, outfile, is_json, override_json=False):
    """
    Download data found in root_url using GET request.
    """
    mode = "wb" if is_image else "w"
    # Get data from the api url
    fullurl = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
    response = requests.get(url=fullurl, params=dict(k=api_key))
    data = response.content if response.status_code == 200 else b""
    # Leave data encoded if its an image
    data = data if is_image else data.decode()
    # Format data and write to outfile.
    if is_json:
        data = get_json_data(json_text=data, override_json=override_json, outfile=outfile)
    open(outfile, mode).write(data)

    return


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


def load_db(uid, stream, pid=None, is_json=False, load_recipe=False, jobs=False, privacy=Project.PRIVATE):
    """
    Load "stream" into database object "uid".
    Loads object as a project by default.
    """
    def project():
        project = Project.objects.get_all(uid=uid).first()
        if not project:
            # Create empty object if not present and populate.
            # Select project owner.
            user = User.objects.filter(is_staff=True).first()
            project = auth.create_project(user=user, name="Project Name", uid=uid, privacy=privacy)
        conf = hjson.loads(stream.read())
        name = conf.get("settings", {}).get("name", project.name)
        text = conf.get("settings", {}).get("help", project.text)
        Project.objects.get_all(uid=uid).update(name=name, text=text)

        return project

    def recipe():
        recipe = Analysis.objects.get_all(uid=uid).first()
        project = Project.objects.get_all(uid=pid).first()
        if not recipe:
            # Create empty object if not present then populate.
            if not project:
                logger.error(f"*** Project id:{pid} does not exist.")
                logger.error(f"\n*** Run `manage.py api project -load --pid={pid}` to create it.")
                sys.exit()
            recipe = auth.create_analysis(project=project, json_text="", template="", uid=uid, name="Recipe Name")

        if is_json:
            data = hjson.loads(stream.read())
            name = data.get("settings", {}).get("name", recipe.name)
            text = data.get("settings", {}).get("help", recipe.text)
            Analysis.objects.get_all(uid=uid).update(json_text=hjson.dumps(data), name=name, text=text)
        else:
            Analysis.objects.get_all(uid=uid).update(template=stream.read())

        if jobs:
            # When creating a job automatically for data in projects
            # it will try to match the value of the parameter to the data name.
            missing_name = ''
            for key, obj in recipe.json_data.items():
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
        return recipe

    return recipe() if load_recipe else project()


def upload(uid, root_dir, pid=None, root_url=None, api_key="", view="recipe_api_template", fname="",
           is_image=False, load_recipe=False, is_json=False, privacy=Project.PRIVATE, jobs=False):

    """
    Upload data into a remote host using API.
    Defaults to local database if root_url is None.
    """

    target = os.path.join(root_dir, uid, fname)
    mode = "rb" if is_image else "r"
    if not os.path.exists(target):
        stream = open(get_thumbnail(), mode) if is_image else io.StringIO("")
    else:
        stream = open(target, mode)
    # Upload to remote host when url is set.
    if root_url:
        return remote_upload(stream=stream, root_url=root_url, uid=uid, api_key=api_key, view=view)
    # Update database info
    if is_image:
        # Update image file.
        mtype = Analysis if load_recipe else Project
        obj = mtype.objects.get_all(uid=uid).first()
        return change_image(obj=obj, file_object=stream)

    return load_db(uid=uid, pid=pid, stream=stream, is_json=is_json, load_recipe=load_recipe, privacy=privacy,
                   jobs=jobs)


def get_data_placeholder(is_json, is_image, uid):

    if is_image:
        placeholder = open(get_thumbnail(), "rb").read()
    elif is_json:
        data = dict(settings=dict(uid=uid,
                                  name="Object Name",
                                  image=f"{uid}.png",
                                  help="Help Text"))
        placeholder = hjson.dumps(data)
    else:
        placeholder = ""

    return placeholder


def download(uid, root_dir, root_url=None, api_key="", is_json=False, view="recipe_api_template",
             fname="", is_image=False, mtype=Analysis, override_json=False):

    # Get placeholder in case object has no image.
    img_path = lambda o: o.image.path if o.image else get_thumbnail()
    mode = "wb" if is_image else "w"
    # Make output directory.
    outdir = os.path.join(root_dir, uid)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, fname)
    if root_url:
        remote_download(root_url=root_url, api_key=api_key, view=view, uid=uid,
                        is_image=is_image, outfile=outfile, is_json=is_json, override_json=override_json)
        return
    # Get data from database
    obj = mtype.objects.get_all(uid=uid).first()
    if not obj:
        data = get_data_placeholder(is_json=is_json, is_image=is_image, uid=uid)
        open(outfile, mode).write(data)
        return
    if is_image:
        data = open(img_path(obj), "rb").read()
    elif is_json:
        data = get_json_data(json_text=obj.json_text, override_json=override_json, outfile=outfile)
    else:
        data = obj.template

    open(outfile, mode).write(data)
    return outfile


def get_recipes(pid, root_url=None, api_key="", rid=""):
    """
    Return recipes belonging to project 'pid' from api if 'root_url' is given
    else return from database.
    """
    # Filter remote site results by 'pid'
    filter_func = lambda key: recipes[key]["project_uid"] == pid
    # Filter by 'rid' instead if that is given.
    if rid:
        filter_func = lambda key: key == rid
    if root_url:
        # Get the recipes from remote url.
        recipe_api = build_api_url(root_url=root_url, api_key=api_key)
        recipes = hjson.loads(requests.get(url=recipe_api, params=dict(k=api_key)).content)
        # Filter recipes from remote host.
        return list(filter(filter_func, recipes))
    query = Q(uid=rid) if rid else Q(project__uid=pid)
    recipes = Analysis.objects.get_all().filter(query)
    if recipes:
        recipes = recipes.values_list("uid", flat=True)
    else:
        # Allows for the creation of 'rid' if it doesn't exist.
        recipes = [rid] if rid else []
    return recipes


def get_image_name(uid, root_url=None, json="conf.hjson", root_dir=None, api_key="", view="recipe_api_json",
                   mtype=Analysis):

    # Get json from url
    if root_url:
        fullurl = build_api_url(root_url=root_url, view=view, uid=uid, api_key=api_key)
        response = requests.get(url=fullurl, params=dict(k=api_key))
        json_text = response.text if response.status_code == 200 else ""
    # Get json from a file
    elif root_dir:
        path = os.path.join(root_dir, uid, json)
        json_text = open(path).read() if os.path.exists(path) else ""
    # Get json from database
    else:
        obj = mtype.objects.get_all(uid=uid).first()
        json_text = obj.json_text if obj else ""

    json_settings = hjson.loads(json_text).get("settings", {})
    # Get the local image name from "settings" in json.
    # Defaults to uid.png
    name = json_settings.get("image", f"{uid}.png")

    return name


def recipe_loader(root_dir, pid, api_key="", root_url=None, rid="", jobs=False):
    """
        Load recipes into api/database from a project found in project_dir.
        Uses PUT request so 'api_key' is required with 'root_url'.
    """
    if not os.path.exists(root_dir):
        logger.error(f"*** Project directory: {root_dir} does not exist.")
        sys.exit()

    # Every subdir in 'project_dir' is a recipe_dir.
    recipe_dirs = [r.name for r in os.scandir(root_dir) if r.is_dir()]
    # Get the specific recipe to load if given.
    recipe_dirs = list(filter(lambda recipe_uid: recipe_uid == rid, recipe_dirs)) if rid else recipe_dirs
    # Prepare the main function used to load.
    load = partial(upload, root_dir=root_dir, root_url=root_url, api_key=api_key, pid=pid, load_recipe=True)
    # Get image name from conf file in directory
    img = lambda uid: get_image_name(uid=uid, root_dir=root_dir)
    for recipe_uid in recipe_dirs:
        load(uid=recipe_uid, fname="conf.hjson", view="recipe_api_json", is_json=True)
        load(uid=recipe_uid, fname=img(uid=recipe_uid), view="recipe_api_image", is_image=True)
        load(uid=recipe_uid, fname="template.sh", jobs=jobs)

        print(f"*** Loaded recipe id: {recipe_uid}")

    return recipe_dirs


def recipe_dumper(root_dir, pid, root_url=None, api_key="", rid=""):
    """
    Dump recipes from the api/database into a target directory
    belonging to single project.
    """
    # Get the recipes from API or database.
    recipes = get_recipes(pid=pid, root_url=root_url, api_key=api_key, rid=rid)
    dump = partial(download, root_url=root_url, root_dir=root_dir, api_key=api_key, override_json=True)
    # Get image name from json on remote host or local database
    img = lambda uid: get_image_name(uid=uid, root_url=root_url, api_key=api_key)
    for recipe_uid in recipes:
        # Dump json, template, and image for a given recipe
        dump(uid=recipe_uid, fname="conf.hjson", is_json=True, view="recipe_api_json")
        dump(uid=recipe_uid, fname=img(uid=recipe_uid), is_image=True, view="recipe_api_image")
        dump(uid=recipe_uid, fname="template.sh")

        print(f"*** Dumped recipe id: {recipe_uid}")
    return recipes


def project_loader(pid, root_dir, root_url=None, api_key="", data=False, data_root="", privacy=""):
    """
    Load project from root_dir into remote host or local database
    """
    pmap = {"private": Project.PRIVATE, "public": Project.PUBLIC}
    privacy = pmap.get(privacy, Project.PRIVATE)

    # Prepare function used to upload
    load = partial(upload, uid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key, privacy=privacy)
    # Get image name from conf file in directory
    img_name = get_image_name(uid=pid, mtype=Project, root_dir=root_dir, json="conf.hjson")

    load(view="project_api_info", fname="conf.hjson")
    load(is_image=True, view="project_api_image", fname=img_name)

    if data:
        json_file = open(os.path.join(root_dir, pid, "conf.hjson"), "r")
        json_data = hjson.load(json_file).get("data", [])
        data_from_json(root=data_root, pid=pid, json_data=json_data)

    print(f"*** Loaded project ({pid}).")
    return


def project_dumper(pid, root_dir, root_url=None, api_key=""):
    """
    Dump project from remote host or local database into root_dir
    """

    # Prepare function used to download info and images
    dump = partial(download, mtype=Project, uid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)
    # Get image name from json on remote host or database
    img_name = get_image_name(uid=pid, mtype=Project, root_url=root_url, view="project_api_info")

    # Dump the project json and image
    dump(fname="conf.hjson", view="project_api_info", is_json=True)
    dump(fname=img_name, view="project_api_image", is_image=True)

    print(f"*** Dumped project {pid}: {root_dir}.")
    return


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

    return


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

    def manage_project(self, **options):

        self.check_error(**options)
        load = options.get("load")
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir") or os.getcwd()
        pid = options.get("pid")
        data = options.get("data")
        privacy = options.get("privacy")
        data_root = options.get("data_root")

        if load:
            project_loader(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key, data=data,
                           data_root=data_root, privacy=privacy)
        else:
            project_dumper(pid=pid, root_dir=root_dir, root_url=root_url, api_key=api_key)
        return

    def manage_recipe(self, **options):

        self.check_error(**options)
        load = options.get("load")
        root_url = options.get("url")
        api_key = options.get("key")
        root_dir = options.get("dir") or os.getcwd()
        uid = options.get("uid")
        pid = options.get("pid")
        create_job = options.get("jobs")

        project_dir = os.path.join(root_dir, pid)
        if load:
            recipes = recipe_loader(root_dir=project_dir, root_url=root_url, api_key=api_key,
                                    rid=uid, pid=pid, jobs=create_job)
        else:
            recipes = recipe_dumper(root_dir=project_dir, root_url=root_url, api_key=api_key, rid=uid, pid=pid)

        msg = f"{len(recipes)} recipes {'loaded' if load else 'dumped'} into project id:{pid}."
        self.stdout.write(self.style.SUCCESS(msg))
        return

    def manage_data(self, **options):

        uid = options.get("uid")
        pid = options.get("pid")
        path = options.get("path")
        name = options.get("name")
        type = options.get("type")
        update_toc = options.get("update_toc")

        if not (pid or uid):
            logger.error("[error] --pid or --uid need to be set.")
            self.run_from_argv(sys.argv + ["--help"])
            sys.exit()

        data_loader(pid=pid, path=path, uid=uid, update_toc=update_toc, name=name, type=type)

        return

    def add_data_commands(self, parser):

        parser.add_argument("--path", type=str, help="Path to data.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to create data in.")
        parser.add_argument("--uid", type=str, help="Data uid to load or update.")
        parser.add_argument('--text', default='', help="A file containing the description of the data")
        parser.add_argument('--name', default='', help="Sets the name of the data")
        parser.add_argument('--type', default='data', help="Sets the type of the data")
        parser.add_argument("--list", action='store_true', help="Show a data list.")
        parser.add_argument("--update_toc", action="store_true", help="Update table of contents for data --uid.")

    def add_recipe_commands(self, parser):
        self.add_api_commands(parser=parser)
        parser.add_argument("--jobs", action='store_true', help="Also creates a queued job for the recipe")
        parser.add_argument('--uid', type=str, default="", help="Recipe uid to load or dump.")
        parser.add_argument('--dir', default='', help="Directory to store/load recipe from.")
        parser.add_argument("--list", action='store_true', help="Show a recipe list.")

    def add_project_commands(self, parser):
        self.add_api_commands(parser=parser)
        parser.add_argument('--privacy', default="private", help="""Privacy of project, only used when creating.""")
        parser.add_argument('--dir', default='', help="Directory to store/load project from.")
        parser.add_argument("--list", action='store_true', help="Show a project list.")
        parser.add_argument("--data", action='store_true', help="Load data found in conf file to local database.")
        parser.add_argument("--data_root", default="", help="Root directory to data found in conf file.")

    def add_api_commands(self, parser):
        """Add default api commands to parser"""
        parser.add_argument('-l', "--load", action="store_true",
                                                   help="""Load to url from a directory.
                                                        Load to database if --url is not set.""")
        parser.add_argument('-d', "--dump", action="store_true",
                                            help="""Dump from a url to directory. 
                                                    Dump from database if --url is not set.""")
        parser.add_argument('--url', default="", help="Site url.")
        parser.add_argument('--key', default='', help="API key. Required to access private projects.")
        parser.add_argument("--pid", type=str, default="", help="Project uid to load from or dump to.")

        return

    def add_job_commands(self, parser):

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

    def add_arguments(self, parser):
        # Load or dump flags

        subparsers = parser.add_subparsers()

        data_parser = subparsers.add_parser("data", help="Data manager.")
        self.add_data_commands(parser=data_parser)

        recipe_parser = subparsers.add_parser("recipe", help="Recipe manager")
        self.add_recipe_commands(parser=recipe_parser)

        project_parser = subparsers.add_parser("project", help="Project manager.")
        self.add_project_commands(parser=project_parser)

        job_parser = subparsers.add_parser("job", help="Job manager.")
        self.add_job_commands(parser=job_parser)

    def check_error(self, **options):
        """
        Check routine errors for parser
        """

        def exit(msg):
            if msg:
                sys.argv.append("--help")
                self.stdout.write(self.style.NOTICE(msg))
                self.run_from_argv(sys.argv)
                sys.exit()
            return

        load = options.get("load")
        dump = options.get("dump")
        root_url = options.get("url")
        api_key = options.get("key")
        pid = options.get("pid")

        if not pid:
            exit(msg="[error] --pid needs to be set.")
        if not (load or dump):
            exit(msg="[error] Set load (-l) or dump (-d) flag.")
        if load and dump:
            exit(msg="[error] Only one flag can be set.")
        if (root_url and load) and not api_key:
            exit(msg="[error] --key is required when loading data to remote site.")

    def handle(self, *args, **options):

        subcommand = sys.argv[2] if len(sys.argv) > 2 else None

        listing = options.get("list")

        if len(sys.argv) == 2:
            self.stdout.write(self.style.NOTICE("Pick a sub-command"))
            self.run_from_argv(sys.argv + ["--help"])
            sys.exit()

        if listing:
            list_obj(mtype=subcommand)
            return

        if subcommand == "project":
            self.manage_project(**options)
            return

        if subcommand == "data":
            self.manage_data(**options)
            return

        if subcommand == "recipe":
            self.manage_recipe(**options)
            return

        if subcommand == "job":
            self.manage_job(**options)
            return
