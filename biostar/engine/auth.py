import difflib
import logging
import uuid, copy
from mimetypes import guess_type

import hjson
from django.conf import settings
from django.contrib import messages
from django.contrib.messages.storage import fallback
from django.db.models import Q
from django.template import Template, Context
from django.template import loader
from django.test import RequestFactory
from django.utils.safestring import mark_safe

from . import models
from . import util
from .const import *
from .models import Data, Analysis, Job, Project, Access

logger = logging.getLogger("engine")


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def fake_request(url="/", data={}, user="", method="POST"):
    "Make a fake request; defaults to POST."

    methods = {"POST": RequestFactory().post, "GET": RequestFactory().get}

    assert method in methods

    request = methods[method](url, data)

    # Mimic messaging system
    request.session = {}
    messages = fallback.FallbackStorage(request=request)
    request._messages = messages

    request.user = user

    return request


def access_denied_message(user, needed_access):
    """
    Generates the access denied message
    """
    tmpl = loader.get_template('widgets/access_denied_message.html')

    # Get the string format of the access.
    needed_access = dict(Access.ACCESS_CHOICES).get(needed_access)

    context = dict(user=user, needed_access=needed_access)
    return tmpl.render(context=context)


def copy_uid(request, instance, board):
    """
    Used to append instance.uid into request.session[board]
    """

    if instance is None:
        messages.error(request, "Object does not exist.")
        return []

    if request.user.is_anonymous:
        messages.error(request, "You need to be logged in.")
        return []

    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})

    board_items = clipboard.get(board, [])
    board_items.append(instance.uid)
    # No duplicates in clipboard

    clipboard[board] = list(set(board_items))

    request.session.update({settings.CLIPBOARD_NAME: clipboard})

    messages.success(request, f"Copied items, there are {len(set(board_items))} in clipboard.")

    return board_items


def authorize_run(user, recipe):
    """
    Returns runnable.
    """

    if user.is_anonymous:
        return False

    if user.is_staff:
        return True

    # Only users with write access can run recipes
    entry = Access.objects.filter(project=recipe.project, user=user, access=Access.WRITE_ACCESS).first()
    if entry:
        return recipe.runnable()

    return False


def generate_script(job):
    """
    Generates a script from a job.
    """
    work_dir = job.path
    json_data = hjson.loads(job.json_text)

    # The base url to the site.
    url_base = f'{settings.PROTOCOL}://{settings.SITE_DOMAIN}{settings.HTTP_PORT}'

    # Extra context added to the script.
    runtime = dict(
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

    # Add the runtime context to the data.
    json_data['runtime'] = runtime
    try:
        # Generate the script.
        template = Template(job.template)
    except Exception as exc:
        template = Template(f"Error loading script template : \n{exc}.")

    context = Context(json_data)
    script = template.render(context)

    return json_data, script


def template_changed(analysis, template):
    """
    Detects a change in the template.
    """
    # Empty template is seen as no change ( False)
    text1 = template.splitlines(keepends=True)
    text2 = analysis.template.splitlines(keepends=True)

    change = list(difflib.unified_diff(text1, text2))

    # print(f"Change: {bool(change)} {change}")
    return change


def get_project_list(user, include_public=True):
    """
    Return projects visible to a user.
    """
    privacy = None
    if include_public:
        privacy = Project.PUBLIC

    if user.is_anonymous:
        # Unauthenticated users see public projects.
        cond = Q(privacy=Project.PUBLIC)
    else:
        # Authenticated users see public projects and private projects with access rights.
        cond = Q(owner=user, privacy=Project.PRIVATE) | Q(privacy=privacy) | Q(access__user=user,
                                                                               access__access__in=[Access.READ_ACCESS,
                                                                                                   Access.WRITE_ACCESS])
    # Generate the query.
    query = Project.objects.filter(cond).distinct()

    return query


def create_project(user, name, uid=None, summary='', text='', stream=None,
                   privacy=Project.PRIVATE, sticky=False, update=False):

    # Set or create the project uid.
    uid = uid or util.get_uuid(8)

    # Attempts to select the project.
    project = Project.objects.filter(uid=uid)

    # If it is not an update request return the project unchanged.
    if project and not update:
        return project.first()

    if project:
        # Update existing project.
        current = project.first()
        text = text or current.text
        name = name or current.name
        project.update(text=text, name=name)
        project = project.first()
        logger.info(f"Updated project: {project.name} uid: {project.uid}")
    else:
        # Create a new project.
        project = Project.objects.create(
            name=name, uid=uid, text=text, owner=user, privacy=privacy, sticky=sticky)
        logger.info(f"Created project: {project.name} uid: {project.uid}")

    if stream:
        # Update the image for the project.
        project.image.save(stream.name, stream, save=True)

    return project


def create_analysis(project, json_text, template, uid=None, user=None, summary='',
                    name='', text='', stream=None, sticky=False, security=Analysis.UNDER_REVIEW, update=False):
    owner = user or project.owner

    analysis = Analysis.objects.get_all(uid=uid)

    # Only update when there is a flag
    if analysis and update:
        # Update analysis
        current = analysis.first()
        text = text or current.text
        name = name or current.name
        template = template or current.template
        json_text = json_text or current.json_text
        analysis.update(text=text, name=name, template=template, json_text=json_text)
        analysis = analysis.first()
        logger.info(f"Updated analysis: uid={analysis.uid} name={analysis.name}")
    else:
        # Create a new analysis
        uid = None if analysis else uid
        text = summary + "\n" + text
        analysis = Analysis.objects.create(project=project, uid=uid, json_text=json_text,
                                           owner=owner, name=name, text=text, security=security,
                                           template=template, sticky=sticky)
        logger.info(f"Created analysis: uid={analysis.uid} name={analysis.name}")

    if stream:
        analysis.image.save(stream.name, stream, save=True)

    return analysis


def make_job_title(recipe, data):
    """
    Creates informative job title that shows job parameters.
    """
    params = data.values()

    # Extracts the field that gets displayed for a parameter
    def extract(param):
        if not param.get("display"):
            return None
        if param.get("source"):
            return param.get("name")
        return param.get("value")

    vals = map(extract, params)
    vals = filter(None, vals)
    vals = map(str, vals)
    vals = ", ".join(vals)

    if vals:
        name = f"Results for {recipe.name}: {vals}"
    else:
        name = f"Results for {recipe.name}"

    return name


def create_job(analysis, user=None, json_text='', json_data={}, name=None, state=Job.QUEUED, uid=None, save=True):
    state = state or Job.QUEUED
    owner = user or analysis.project.owner
    project = analysis.project

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or analysis.json_text

    # Needs the json_data to set the summary.
    json_data = hjson.loads(json_text)

    # Generate a meaningful job title.
    name = make_job_title(recipe=analysis, data=json_data)

    # Create the job instance.
    job = Job(name=name, state=state, json_text=json_text,
              security=Job.AUTHORIZED, project=project, analysis=analysis, owner=owner,
              template=analysis.template, uid=uid)

    if save:
        job.save()
        logger.info(f"Created job id={job.id} name={job.name}")

    return job


def delete_object(obj, request):
    obj.deleted = not obj.deleted
    obj.save()
    msg = f"Deleted <b>{obj.name}</b>." if obj.deleted else f"Restored <b>{obj.name}</b>."
    messages.success(request, mark_safe(msg))

    return obj.deleted


def has_write_access(user, project):
    """
    Returns True if a user has write access to an instance
    """

    # Anonymous user may not have write access.
    if user.is_anonymous:
        return False

    # Users that may access a project.
    cond1 = (user == project.owner) or user.is_staff

    # User has been given write access to the project
    cond2 = models.Access.objects.filter(user=user, project=project,
                                         access=models.Access.WRITE_ACCESS).first()

    # One of the conditions has to be true.
    access = cond1 or cond2

    return access


def guess_mimetype(fname):
    "Return mimetype for a known text filename"

    mimetype, encoding = guess_type(fname)

    ext = os.path.splitext(fname)[1].lower()

    # Known text extensions ( .fasta, .fastq, etc.. )
    if ext in KNOWN_TEXT_EXTENSIONS:
        mimetype = 'text/plain'

    return mimetype


def create_path(fname, data):
    """
    Returns a proposed path based on fname to the storage folder of the data.
    Attempts to preserve the extension but also removes all whitespace from the filenames.
    """
    # Select the file name.

    fname = os.path.basename(fname)

    # The data storage directory.
    data_dir = data.get_data_dir()

    # Make the data directory if it does not exist.
    os.makedirs(data_dir, exist_ok=True)

    # Build the file name under the new location.
    path = os.path.abspath(os.path.join(data_dir, fname))

    return path


def link_file(path, data):
    dest = create_path(fname=path, data=data)

    if not os.path.exists(dest):
        os.symlink(path, dest)

    return dest


def fill_data_by_name(project, json_data):
    """
    Fills json information by name. Used when filling in
    demonstration data and not user selection.
    """

    json_data = copy.deepcopy(json_data)
    # A mapping of data by name
    store = dict((data.name, data) for data in project.data_set.all())

    for field, item in json_data.items():
        # If the field is a data field then fill in more information.
        if item.get("source") == "PROJECT":
            name = item.get("value")
            data = store.get(name)
            if data:
                # This mutates the `item` dictionary!
                data.fill_dict(item)

    return json_data


def create_data(project, user=None, stream=None, path='', name='',
                text='', type='', uid=None, summary=''):
    # We need absolute paths with no trailing slashes.
    path = os.path.abspath(path).rstrip("/")

    # Create the data.
    type = type or "DATA"

    # Set the object unique id.
    uid = uid or util.get_uuid(8)

    # The owner of the data will be the first admin user if not set otherwise.
    owner = user or models.User.objects.filter(is_staff=True).first()

    # Create the data object.
    data = Data.objects.create(name=name, owner=owner, state=Data.PENDING, project=project,
                               type=type, text=text, uid=uid)

    # The source of the data is a stream and will become the path
    # that should be added.
    if stream:
        name = name or stream.name
        fname = '_'.join(name.split())
        path = create_path(data=data, fname=fname)
        util.write_stream(stream=stream, dest=path)
        # Mark incoming file as uploaded
        data.method = Data.UPLOAD

    # The path is a file.
    isfile = os.path.isfile(path)

    # The path is a directory.
    isdir = os.path.isdir(path)

    # The data is a single file on a path.
    if isfile:
        link_file(path=path, data=data)
        logger.info(f"Linked file: {path}")

    # The data is a directory.
    # Link each file of the directory into the storage directory.
    if isdir:
        for p in os.scandir(path):
            link_file(path=p.path, data=data)
            logger.info(f"Linked file: {p}")

    # Invalid paths and empty streams still create the data
    # but set the data state will be set to error.
    missing = not(path or stream)

    # An invalid entry here.
    if path and missing:
        state = Data.ERROR
        logger.error(f"Invalid data path: {path}")
    else:
        state = Data.READY

    # Make the table of content file with files in data_dir.
    data.make_toc()

    # Set updated attributes
    data.state = state
    data.name = name or os.path.basename(path) or 'Data'

    # Trigger another save.
    data.save()

    # Set log for data creation.
    logger.info(f"Added data type={data.type} name={data.name} pk={data.pk}")

    return data
