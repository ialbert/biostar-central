import difflib
import logging
import uuid
from mimetypes import guess_type

import hjson
from biostar import settings
from biostar.accounts.models import Profile
from django.contrib import messages
from django.db.models import Q
from django.template import Template, Context
from django.template import loader
from django.test.client import RequestFactory
from django.utils.safestring import mark_safe

from . import util
from .const import *
from .models import Data, Analysis, Job, Project, Access


logger = logging.getLogger("engine")


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


def join(*args):
    return os.path.abspath(os.path.join(*args))


def access_denied_message(user, access):
    """
    Generates the access denied message
    """
    tmpl = loader.get_template('widgets/access_denied_message.html')
    context = dict(user=user, access=access)
    return tmpl.render(context=context)


def authorize_analysis(user, recipe):

    if user.is_staff:
        return Analysis.AUTHORIZED

    return Analysis.UNDER_REVIEW

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
        cond = Q(privacy=privacy) | Q(access__user=user, access__access__gt=Access.NO_ACCESS)

    # Generate the query.
    query = Project.objects.filter(cond).distinct()

    return query


def check_obj_access(user, instance, access=Access.NO_ACCESS, request=None, login_required=False, role=None):
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

    project = instance.project

    # Check for logged in user and login requirement.
    if user.is_anonymous and login_required:
        messages.error(request, "You must be logged in to perform that action.")
        return False

    # If the project is public then any user can have read access.
    if project.privacy == Project.PUBLIC and access == Access.READ_ACCESS:
        return True

    # Check ownership access.
    if access == Access.OWNER_ACCESS:
        if project.owner == user or instance.owner == user:
            return True
        else:
            msg = "Only the creator of the object or project can perform that action."
            messages.error(request, msg)
            return False

    # Prepare the access denied message.
    access_text = Access.ACCESS_MAP.get(access, 'Invalid')
    deny = access_denied_message(user=user, access=access_text)

    # Anonymous users have no other access permissions.
    if user.is_anonymous:
        messages.error(request, deny)
        return False

    # Check user access.
    entry = Access.objects.filter(user=user, project=project).first()
    entry = entry or Access(access=Access.NO_ACCESS)
    # Check user role.
    has_role = (role is not None) and user.profile.role == role

    if entry.access < access and not has_role:
        messages.warning(request, deny)
        return False
    if entry.access >= access or has_role:
        return True

    messages.error(request, "Invalid fall through.")
    return False


def create_project(user, name, uid=None, summary='', text='', stream=None,
                   privacy=Project.PRIVATE, sticky=True, update=False):
    uid = uid or util.get_uuid(8)
    project = Project.objects.filter(uid=uid)

    if project and not update:
        return project.first()

    if project:
        # Update project
        current = project.first()
        summary = summary or current.summary
        text = text or current.text
        name = name or current.name
        project.update(summary=summary, text=text, name=name)
        project = project.first()
        logger.info(f"Updated project: {project.name} uid: {project.uid}")
    else:
        project = Project.objects.create(
            name=name, uid=uid, summary=summary, text=text, owner=user, privacy=privacy, sticky=sticky)
        logger.info(f"Created project: {project.name} uid: {project.uid}")

    if stream:
        project.image.save(stream.name, stream, save=True)

    return project


def create_analysis(project, json_text, template, uid=None, user=None, summary='',
                    name='', text='', stream=None, sticky=False, security=Analysis.UNDER_REVIEW, update=False):
    owner = user or project.owner

    analysis = Analysis.objects.filter(uid=uid)

    # Only update when there is a flag
    if analysis and update:
        # Update analysis
        current = analysis.first()
        summary = summary or current.summary
        text = text or current.text
        name = name or current.name
        template = template or current.template
        json_text = json_text or current.json_text
        analysis.update(summary=summary, text=text, name=name, template=template, json_text=json_text)
        analysis = analysis.first()
        logger.info(f"Updated analysis: uid={analysis.uid} name={analysis.name}")
    else:
        # Create a new analysis
        uid = None if analysis else uid
        analysis = Analysis.objects.create(project=project, uid=uid, summary=summary, json_text=json_text,
                                           owner=owner, name=name, text=text, security=security,
                                           template=template, sticky=sticky)
        logger.info(f"Created analysis: uid={analysis.uid} name={analysis.name}")

    if stream:
        analysis.image.save(stream.name, stream, save=True)

    return analysis


def make_job_summary(data, summary='', title='', name="widgets/job_summary.html"):
    '''
    Summarizes job parameters.
    '''

    context = dict(data=data, summary=summary, title=title)
    template = loader.get_template(name)
    result = template.render(context)

    return result


def make_job_title(recipe, data):
    name = f"Results: {recipe.name}"
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
    summary = make_job_summary(json_data, summary=analysis.summary)

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


def check_data_name(name, data, bool=False):

    copy = name
    i = 0
    check_name = lambda n: Data.objects.exclude(pk=data.pk).filter(name=n,
                                                  project=data.project).exists()
    if bool:
        return check_name(name)

    while check_name(copy) and i < 1000:
        i += 1
        copy = f"{name} ({i})"

    return  copy


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

    os.symlink(path, dest)
    return dest


def create_data(project, user=None, stream=None, path='', name='',
                text='', summary='', type="",  uid=None, paths=[]):

    # We need absolute paths with no trailing slashes.
    path = os.path.abspath(path).rstrip("/") if path else ""

    # Create the data.
    type = type or "DATA"
    uid = uid or util.get_uuid(8)
    data = Data.objects.create(name=name, owner=user, state=Data.PENDING, project=project,
                               type=type, summary=summary, text=text, uid=uid)

    # The source of the data is a stream is written into the destination.
    if stream:
        name = name or stream.name
        fname = '_'.join(name.split())
        dest = create_path(data=data, fname=fname)
        util.write_stream(stream=stream, dest=dest)
        # Mark incoming file as uploaded
        data.method = Data.UPLOAD

    # Link single file
    if path:
        link_file(path=path, data=data)
        logger.info(f"Linked file: {path}")

    # Link list of files
    for p in paths:
        link_file(path=p, data=data)
        logger.info(f"Linked file: {p}")

    # Invalid paths and empty streams still create the data but set the data state to error.
    missing = not (os.path.isdir(path) or os.path.isfile(path) or stream)
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
    data.summary = summary

    # Trigger another save.
    data.save()

    # Set log for data creation.
    logger.info(f"Added data type={data.type} name={data.name} pk={data.pk}")

    return data
