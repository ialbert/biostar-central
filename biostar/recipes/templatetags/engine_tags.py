
import logging
import os
import re
from textwrap import dedent
import hashlib
import urllib.parse
import random

from datetime import timedelta, datetime
from django.contrib import messages
from django import template
from django.conf import settings
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.core.paginator import Paginator
from django.db.models import Q,Count
from django.template import defaultfilters
from django.utils.safestring import mark_safe

from biostar.recipes import auth, util, const
from biostar.recipes.models import Job, make_html, Project, Data, Analysis, Access, CommandType, Command


logger = logging.getLogger("engine")
register = template.Library()

JOB_COLORS = {
    Job.SPOOLED: "violet",
    Job.ERROR: "red", Job.QUEUED: "teal",
    Job.RUNNING: "orange", Job.COMPLETED: "green"
}

DATA_COLORS = {
    Data.PENDING: "teal", Data.READY: "green", Data.ERROR: "red"
}

@register.simple_tag
def randparam():
    "Append to URL to bypass server caching of CSS or JS files"
    return f"?randval={random.randint(1, 10000000)}" if settings.DEBUG else ""

@register.filter
def time_ago(date):
    pluralize = lambda value, word: f"{value} {word}s" if value > 1 else f'{value} {word}'
    if not date:
        return ''
    delta = util.now() - date
    if delta < timedelta(minutes=1):
        return 'just now'
    elif delta < timedelta(hours=1):
        unit = pluralize(delta.seconds // 60, "minute")
    elif delta < timedelta(days=1):
        unit = pluralize(delta.seconds // 3600, "hour")
    elif delta < timedelta(days=30):
        unit = pluralize(delta.days, "day")
    elif delta < timedelta(days=90):
        unit = pluralize(int(delta.days / 7), "week")
    elif delta < timedelta(days=730):
        unit = pluralize(int(delta.days / 30), "month")
    else:
        diff = delta.days / 365.0
        unit = '%0.1f years' % diff
    return "%s ago" % unit

def join(*args):
    return os.path.abspath(os.path.join(*args))

@register.filter
def bignum(number):
    "Reformats numbers with qualifiers as K"
    try:
        value = float(number) / 1000.0
        if value > 10:
            return "%0.fk" % value
        elif value > 1:
            return "%0.1fk" % value
    except ValueError as exc:
        pass
    return str(number)


@register.simple_tag()
def user_score(user):
    score = user.profile.score * 10
    return score


@register.inclusion_tag('widgets/user_icon.html')
def user_icon(user):
    score = user_score(user)
    context = dict(user=user, score=score)
    return context


@register.inclusion_tag('widgets/list_view.html', takes_context=True)
def list_projects(context, target):

    user = context["request"].user
    projects = auth.get_project_list(user=target)

    # Don't show private projects non owners
    if user != target:
        projects = projects.exclude(privacy=Project.PRIVATE)

    projects = projects.annotate(data_count=Count('data', distinct=True, filter=Q(deleted=False)),
                                 job_count=Count('job', distinct=True, filter=Q(deleted=False)),
                                 recipe_count=Count('analysis', distinct=True, filter=Q(deleted=False)),
                                 )
    projects = projects.order_by("-rank", "-lastedit_date")

    return dict(projects=projects)


@register.simple_tag
def gravatar(user, size=80):
    style = "retro"
    if user.is_anonymous or user.profile.is_suspended:
        # Removes spammy images for suspended users
        # email = 'suspended@biostars.org'.encode('utf8')
        style = "monsterid"
    else:
        if user.profile.is_moderator:
            style = "robohash"
        email = user.email.encode('utf8')

    hash = hashlib.md5(email).hexdigest()

    gravatar_url = "https://secure.gravatar.com/avatar/%s?" % hash
    gravatar_url += urllib.parse.urlencode({
        's': str(size),
        'd': style,
    }
    )
    return gravatar_url


@register.filter
def highlight(source, target):

    # Look for case insensitive matches in the source
    highlighting = re.search(f"(?i){target}", source)

    target = highlighting.group() if highlighting else target

    # Highlight the target
    highlighted = mark_safe(f"<div class='match'>{target}</div>")

    return source.replace(target, highlighted)


@register.simple_tag
def get_qiime2view_link(file_serve_url):
    port = f':{settings.HTTP_PORT}' if settings.HTTP_PORT else ''
    site = f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}{port}"

    full_url = site + file_serve_url

    qiime_link = util.qiime2view_link(full_url)

    return qiime_link


@register.inclusion_tag('widgets/list_view.html', takes_context=True)
def list_view(context, projects=None, data_list=None, recipe_list=None, job_list=None):
    request = context["request"]

    return dict(projects=projects, data_list=data_list, recipe_list=recipe_list,
                job_list=job_list, request=request)


@register.inclusion_tag('widgets/paste.html', takes_context=True)
def paste(context, project, current=""):

    request = context["request"]
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})
    board = clipboard.get(current, [])

    clipboard_count = len(board) if request.user.is_authenticated else 0

    extra_context = dict(clipboard_count=clipboard_count, project=project, current=current, board=board, context=context)
    context.update(extra_context)
    return context


@register.filter
def is_checkbox(field):
    "Check if current field is a checkbox"

    try:
        if field.field.widget.input_type == "checkbox":
            return True
    except Exception as exc:
        logger.error(exc)

    return False


@register.filter
def is_qiime_archive(file=None):
    filename = file if isinstance(file, str) else file.path

    return filename.endswith(".qza") or filename.endswith(".qzv")


@register.simple_tag
def privacy_label(project):
    label = mark_safe(f'<span class ="ui label">{project.get_privacy_display()}</span>')
    return label


@register.inclusion_tag('widgets/authorization_required.html', takes_context=True)
def security_label(context, analysis):
    context.update(dict(analysis=analysis))

    return context


@register.simple_tag
def job_color(job):
    """
    Returns a color based on job status.
    """
    return JOB_COLORS.get(job.state, "")


@register.simple_tag
def activate(value1, value2):
    """
    Returns a color based on job status.
    """
    return "active" if value1 == value2 else ''


@register.simple_tag
def data_color(data):
    "Return a color based on data status."

    return DATA_COLORS.get(data.state, "")


@register.simple_tag
def type_label(data):
    if data.type:
        label = lambda x: f"<span class='ui label' > {x} </span>"
        types = [label(t) for t in data.type.split(',')]
        return mark_safe(''.join(types))
    return ""


@register.simple_tag
def state_label(data, error_only=False):
    label = f'<span class="ui { DATA_COLORS.get(data.state, "") } label"> {data.get_state_display()} </span>'

    # Error produce error only.
    if error_only and data.state not in (Data.ERROR, Data.PENDING):
        label = ""

    return mark_safe(label)


@register.simple_tag
def img(obj):
    """
    Returns the image associated with the object or a placeholder
    """
    if obj.image:
        return obj.image.url
    else:
        return static("images/placeholder.png")


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.inclusion_tag('widgets/project_title.html', takes_context=True)
def project_title(context, project):
    """
    Returns a label for project.
    """
    return dict(project=project)


@register.inclusion_tag('widgets/recipe_form.html')
def recipe_form(form):
    """
    Renders a recipe form.
    """
    return dict(form=form)


@register.inclusion_tag('widgets/interface_options.html')
def interface_options():
    return dict()


@register.inclusion_tag('widgets/code_options.html', takes_context=True)
def code_options(context):

    command_types = CommandType.objects.all()
    extra_context = dict(command_types=command_types)
    context.update(extra_context)
    return context


@register.inclusion_tag('widgets/json_field.html')
def json_field(json_text):

    context = dict(json_text=json_text)
    return context


@register.inclusion_tag('widgets/template_field.html')
def template_field(tmpl):

    context = dict(template=tmpl)
    return context


@register.inclusion_tag('widgets/created_by.html')
def created_by(date, user=None):
    """
    Renders a created by link
    """
    return dict(date=date, user=user)

@register.inclusion_tag('widgets/access_form.html')
def access_form(project, user, form):
    """
    Generates an access form.
    """

    return dict(project=project, user=user, form=form)


@register.inclusion_tag('widgets/job_elapsed.html')
def job_minutes(job):

    check_back = 'check_back' if job.state in [Job.SPOOLED, Job.RUNNING] else ''

    return dict(job=job, check_back=check_back)


@register.simple_tag
def size_label(data):
    """
    Returns a label for data sizes.
    """

    size = f"{defaultfilters.filesizeformat(data.size)}"
    return mark_safe(f"<span class='ui mini label'>{size}</span>")


@register.inclusion_tag('widgets/directory_list.html', takes_context=True)
def directory_list(context, obj):
    """
    Generates an HTML listing for files in a directory.
    """

    # Starting location.
    root = obj.get_data_dir()

    # The serve url depends on data type..
    serve_url = "job_serve" if isinstance(obj, Job) else "data_serve"
    copy_url = "job_file_copy" if isinstance(obj, Job) else "data_file_copy"

    # This will collet the valid filepaths.
    paths = []
    try:
        # Walk the filesystem and collect all files.
        for fpath, fdirs, fnames in os.walk(root, followlinks=True):
            paths.extend([join(fpath, fname) for fname in fnames])
        # Image extension types.
        IMAGE_EXT = {"png", "jpg", "gif", "jpeg"}

        # Add more metadata to each path.
        def transform(path):
            tstamp = os.stat(path).st_mtime
            size = os.stat(path).st_size
            rel_path = os.path.relpath(path, root)
            elems = os.path.split(rel_path)
            dir_names = elems[:-1]
            if dir_names[0] == '':
                dir_names = []
            last_name = elems[-1]
            is_image = last_name.split(".")[-1] in IMAGE_EXT
            return rel_path, dir_names, last_name, tstamp, size, is_image

        # Transform the paths.
        paths = map(transform, paths)

        # Sort by the tuple fields..
        paths = sorted(paths)

    except Exception as exc:
        logging.error(exc)
        paths = []

    return dict(paths=paths, obj=obj, serve_url=serve_url, copy_url=copy_url, user=context["request"].user)


@register.inclusion_tag('widgets/form_errors.html')
def form_errors(form):
    """
    Turns form errors into a data structure
    """
    try:
        errorlist = [('', message) for message in form.non_field_errors()]

        for field in form:
            for error in field.errors:
                errorlist.append((f'{field.name}:', error))
    except Exception:
        errorlist = []

    context = dict(errorlist=errorlist)

    return context


@register.filter
def markdown(text):
    """
    Generates HTML from a markdown value.
    """
    if not text:
        return ''
    text = dedent(text)
    html = make_html(text)
    return mark_safe(html)


@register.inclusion_tag("widgets/menubar.html", takes_context=True)
def menubar(context, request=None, with_search=True):
    user = context.request.user
    context.update(dict(user=user, request=request, with_search=with_search))

    return context
