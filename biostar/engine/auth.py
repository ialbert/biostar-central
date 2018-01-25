import difflib
import logging
import uuid

import hjson
from django.contrib import messages
from django.db.models import Q
from django.template import Template, Context
from django.template import loader
from django.test.client import RequestFactory
from django.utils.safestring import mark_safe

from biostar.engine import settings
from . import util
from .const import *
from .models import Data, Analysis, Job, Project, Access
from .templatetags import engine_tags

CHUNK = 1024 * 1024

logger = logging.getLogger("engine")


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def get_analysis_attr(analysis, project=None):
    "Get analysis attributes"

    project = project or analysis.project
    json_text = analysis.json_text
    template = analysis.template
    owner = analysis.owner
    summary = analysis.summary
    name = analysis.name
    text = analysis.text

    return dict(project=project, json_text=json_text, template=template,
                user=owner, summary=summary, name=name, text=text)


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

    # print(f"Change: {bool(change)}")
    return change


def get_project_list(user):
    """
    Return projects visible to a user.
    """

    if user.is_anonymous:
        # Unauthenticated users see public projects.
        cond = Q(privacy=Project.PUBLIC)
    else:
        # Authenticated users see public projects and private projects with access rights.
        cond =  Q(privacy=Project.PUBLIC) | Q(access__user=user, access__access__gt=Access.NO_ACCESS)

    # Generate the query.
    query = Project.objects.filter(cond).distinct()

    return query


def check_obj_access(user, instance, access=Access.ADMIN_ACCESS, request=None, login_required=False):
    """
    Validates object access.
    """
    # This is so that we can inform users in more granularity
    # but also to allow this function to be called outside a web view
    # where there are no requests.
    request = request or RequestFactory()

    # The object does not exist.
    if not instance:
        messages.error(request, "Object does not exist.")
        return False

    # A textual representation of the access
    access_text = Access.ACCESS_MAP.get(access, 'Invalid')

    # Works for projects or objects with an attribute of project.
    if hasattr(instance, "project"):
        project = instance.project
    else:
        project = instance

    # Check for logged in user and login requirement.
    if user.is_anonymous() and login_required:
        msg = f"""
            You must be logged in to perform that action.
        """
        msg = mark_safe(msg)
        messages.error(request, msg)
        return False

    # If the project is public or shareable project and a user is asking for read access.
    if (project.privacy in (Project.PUBLIC, Project.SHAREABLE)):
        if (access in (Access.READ_ACCESS, Access.RECIPE_ACCESS)):
            return True

    # Prepare the access denied message.
    deny = engine_tags.access_denied_message(user=user, access=access_text)

    # Anonymous users have no other access permissions.
    if user.is_anonymous():
        messages.error(request, deny)
        return False

    # Check user access.
    entry = Access.objects.filter(user=user, project=project).first()

    # No access permissions for the user on the project.
    if not entry:
        messages.error(request, deny)
        return False

    # The stored access is less than the required access.
    if entry.access < access:
        messages.warning(request, deny)
        return False

    # Permissions granted to the object.
    if entry.access >= access:
        return True

    # This should never trigger and is here to catch bugs.
    messages.error(request, "Access denied! Invalid fall-through!")
    return False


def get_data(user, project, query, data_type=None):
    """
    Returns a dictionary keyed by data stored in the project.
    """
    if data_type:
        query = query.filter(type=data_type)
    datamap = dict((obj.id, obj) for obj in query)
    logger.info(f"{user.email} got {len(datamap)} data from {project.name}")

    return datamap


def create_project(user, name, uid=None, summary='', text='', stream='',
                   privacy=Project.PRIVATE, sticky=True):
    project = Project.objects.create(
        name=name, uid=uid, summary=summary, text=text, owner=user, privacy=privacy, sticky=sticky)

    if stream:
        project.image.save(stream.name, stream, save=True)

    logger.info(f"Created project: {project.name} uid: {project.uid}")

    return project


def create_analysis(project, json_text, template, json_file=None, uid=None, user=None, summary='',
                    name='', text='', stream=None, sticky=False, security=Analysis.UNDER_REVIEW):
    owner = user or project.owner
    name = name or 'Analysis name'
    text = text or 'Analysis text'

    analysis = Analysis.objects.create(project=project, uid=uid, summary=summary, json_text=json_text,
                                       owner=owner, name=name, text=text, security=security,
                                       template=template, sticky=sticky, json_file=json_file)
    if stream:
        analysis.image.save(stream.name, stream, save=True)

    logger.info(f"Created analysis: uid={analysis.uid} name={analysis.name}")

    return analysis


def make_summary(data, summary='', title='', name="widgets/job_summary.html"):
    '''
    Summarizes job parameters.
    '''

    context = dict(data=data, summary=summary, title=title)
    template = loader.get_template(name)
    result = template.render(context)

    return result


def make_job_title(recipe, data):

    collect = []
    for key, obj in data.items():
        if obj.get('label'):
            value = obj.get('name') or obj.get('value')
            collect.append(str(value))

    if collect:
        label = ", ".join(collect)
        name = f"{recipe.name} on {label}"
    else:
        name = f"{recipe.name}"

    return name


def create_job(analysis, user=None, json_text='', json_data={}, name=None, state=None, uid=None, save=True):
    state = state or Job.QUEUED
    owner = user or analysis.project.owner

    project = analysis.project

    if json_data:
        json_text = hjson.dumps(json_data)
    else:
        json_text = json_text or analysis.json_text

    # Needs the json_data to set the summary.
    json_data = hjson.loads(json_text)

    # Generate the summary from the data.
    summary = make_summary(json_data, summary=analysis.summary)

    # Generate a meaningful job title.
    name = make_job_title(recipe=analysis, data=json_data)

    # Create the job instance.
    job = Job(name=name, summary=summary, state=state, json_text=json_text,
              security=analysis.security, project=project, analysis=analysis, owner=owner,
              template=analysis.template, uid=uid)

    if save:
        job.save()
        logger.info(f"Created job id={job.id} name={job.name}")

    return job


def findfiles(location, collect, skip=""):
    """
    Returns a list of all files in a directory.
    """

    for item in os.scandir(location):

        if os.path.abspath(item.path) == skip:
            continue

        if item.is_dir():
            findfiles(item.path, collect=collect, skip=skip)
        else:
            collect.append(os.path.abspath(item.path))

    return collect


def link_files(path, data, skip='', summary=''):
    "Param : data (instance) - An instance of a Data object"

    files = []
    if path and os.path.isdir(path):
        # If the path is a directory, symlink all files in the directory.
        logger.info(f"Linking path: {path}")
        collect = findfiles(path, collect=[], skip=skip)
        for fname in os.listdir(path):
            if join(path, fname) == skip:
                continue
            dest = create_path(fname=fname, data=data)
            src = os.path.join(path, fname)
            os.symlink(src, dest)
        # Append directory info to data summary
        data.summary = f'Contains {len(collect)} files. {summary}'
        logger.info(f"Linked {len(collect)} files.")
        files.extend(collect)

    elif path and os.path.isfile(path):
        # The path is a file.
        dest = create_path(path, data=data)
        os.symlink(path, dest)
        logger.info(f"Linked file: {path}")
        files.append(dest)

    return files


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


def create_data(project, user=None, stream=None, path='', name='',
                text='', summary='', type="", skip=""):
    "Param : skip (str) - full file path found in 'path' that will be ignored when linking."

    # We need absolute paths with no trailing slashes.
    path = os.path.abspath(path).rstrip("/") if path else ""

    # Create the data.
    type = type or "DATA"
    #uid = util.get_uuid(8)

    data = Data.objects.create(name=name, owner=user, state=Data.PENDING, project=project,
                               type=type, summary=summary, text=text)#, uid=uid)

    # The source of the data is a stream is written into the destination.
    if stream:
        dest = create_path(data=data, fname=stream.name)
        with open(dest, 'wb') as fp:
            chunk = stream.read(CHUNK)
            while chunk:
                fp.write(chunk)
                chunk = stream.read(CHUNK)

    link_files(path=path, skip=skip, data=data, summary=summary)

    # Invalid paths and empty streams still create the data but set the data state to error.
    missing = not (os.path.isdir(path) or os.path.isfile(path) or stream)
    if path and missing:
        state = Data.ERROR
        logger.error(f"Invalid data path: {path}")
    else:
        state = Data.READY

    # Find all files in the data directory, skipping the table of contents.
    tocname = data.get_path()
    collect = findfiles(data.get_data_dir(), collect=[], skip=tocname)

    # Write the table of contents.
    with open(tocname, 'w') as fp:
        fp.write("\n".join(collect))

    # Find the cumulative size of the files.
    size = 0
    for elem in collect:
        if os.path.isfile(elem):
            size += os.stat(elem, follow_symlinks=True).st_size

    # Finalize the data name
    name = name or os.path.split(path)[1] or 'Data'

    # Set updated attributes
    data.size = size
    data.state = state
    data.name = name
    data.summary = summary
    data.file = tocname

    # Trigger another save.
    data.save()

    # Set log for data creation.
    logger.info(f"Added data type={data.type} name={data.name}")

    return data
