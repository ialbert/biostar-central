import logging
import os
import re
from textwrap import dedent
import hashlib
import urllib.parse
import random
from itertools import count, islice

from datetime import timedelta, datetime
from django.contrib import messages
from django.forms.widgets import CheckboxInput
from django import template, forms
from django.shortcuts import reverse
from django.conf import settings
from django.core.paginator import Paginator
from django.db.models import Q, Count
from django.template import defaultfilters
from django.utils.safestring import mark_safe
from django.utils.timezone import utc

from biostar.recipes import auth, util, const
from biostar.recipes.models import Job, make_html, Project, Data, Analysis, Access, SnippetType, Snippet

logger = logging.getLogger("engine")
register = template.Library()

DATA_COLORS = {
    Data.PENDING: "teal", Data.READY: "green", Data.ERROR: "red"
}


@register.simple_tag
def randparam():
    "Append to URL to bypass server caching of CSS or JS files"
    return f"?randval={random.randint(1, 10000000)}" if settings.DEBUG else ""


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

@register.inclusion_tag('widgets/privacy_label.html')
def privacy_label(project):
    context = dict(project=project)
    return context



@register.inclusion_tag('widgets/list_view.html', takes_context=True)
def list_projects(context, target):
    """List projects belonging to a specific user
    """
    user = context["request"].user
    projects = auth.get_project_list(user=target).filter(owner=target)

    # Don't show private projects non owners
    if user != target:
        projects = projects.exclude(privacy=Project.PRIVATE)

    projects = projects.order_by("-rank", "-lastedit_date")

    return dict(projects=projects, user=user, target=target)


@register.filter
def is_job(obj):
    return isinstance(obj, Job)


@register.simple_tag
def render_script(recipe, tmpl, user):

    # Render the script when the
    if user.is_anonymous:
        return auth.render_script(recipe=recipe, tmpl=tmpl)

    return tmpl


@register.inclusion_tag('widgets/pages.html', takes_context=True)
def pages(context, objs, show_step=True):
    request = context["request"]
    url = request.path
    return dict(objs=objs, url=url, show_step=show_step, request=request)


@register.simple_tag
def access_class(user, project):
    """
    CSS class returned based on access to a project
    """

    if user.is_anonymous:
        return ""

    if user == project.owner:
        return "write_access"

    obj = Access.objects.filter(user=user, project=project).first()

    return obj.access if obj else ""


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
def endswith(string, suffix):

    return string.endswith(suffix)


def find_fragments(source, target, nfrags=3, offset=25):

    # Look for case insensitive matches of target in the source
    matches = re.finditer(f"(?i){target}", source)
    matches = islice(zip(count(1), matches), nfrags)
    fragments = []

    # Collect matches as fragments
    for idx, match in matches:
        # Get left side, right side, and center text of match
        left = match.start(0) - offset
        right = match.end(0) + offset
        text = match.group()

        if left < 0:
            left = 0
        if right > len(source):
            right = len(source)

        fragments.append((left,  right, text))

    return fragments


@register.inclusion_tag('widgets/clipboard.html', takes_context=True)
def clipboard(context, project_uid):

    request = context['request']
    user = request.user
    project = Project.objects.filter(uid=project_uid).first()
    board = auth.recent_clipboard(request=request)
    key, vals = board
    board_count = len(vals)
    movable = key in [const.COPIED_RECIPES, const.COPIED_DATA]
    if project and auth.is_readable(user=user, obj=project) and board_count:
        # Load items into clipboard
        context = dict(count=board_count, board=key, is_recipe=key == const.COPIED_RECIPES,
                       movable=movable)
    else:
        context = dict()

    return context


@register.filter
def highlight(source, target):

    # Number of fragments to show
    nfrags = 2

    # Character offset used to pad highlighted items
    offset = 30

    # Applies the highlighter class to each fragment
    def highlighter(parent, sub):
        return parent.replace(sub, mark_safe(f"<div class='match'>{sub}</div>"))

    # Gather the fragments.
    fragments = find_fragments(source=source, target=target, nfrags=nfrags, offset=offset)

    if fragments:
        result = [highlighter(source[start:end], txt) for start, end, txt in fragments]
        result = "...".join(result)
    else:
        result = source[:offset * 4]

    result += "..." if len(source) > len(result) else ""

    return result


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
    user = request.user
    return dict(projects=projects, user=user, data_list=data_list, recipe_list=recipe_list,
                job_list=job_list, request=request)


@register.filter
def is_checkbox(field):
    "Check if current field is a checkbox"

    try:
        state = isinstance(field.field.widget, CheckboxInput)
        return state
    except Exception as exc:
        logger.error(exc)

    return False


@register.filter
def is_qiime_archive(file=None):
    filename = file if isinstance(file, str) else file.path

    return filename.endswith(".qza") or filename.endswith(".qzv")


@register.inclusion_tag('widgets/authorization_required.html', takes_context=True)
def security_label(context, analysis):
    user = context['request'].user

    if user.is_anonymous:
        is_readable = False
    else:
        is_readable = auth.is_readable(user=user, obj=analysis.project)

    context.update(dict(recipe=analysis, is_readable=is_readable))

    return context


@register.simple_tag
def full_url():
    if settings.HTTP_PORT:
        return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}:{settings.HTTP_PORT}"
    else:
        return f"{settings.PROTOCOL}://{settings.SITE_DOMAIN}"


@register.simple_tag
def job_color(job):
    """
    Returns a color based on job status.
    """
    return auth.job_color(job)


@register.simple_tag
def activate(value1, value2):
    """
    Returns a color based on job status.
    """
    return "orange active" if value1 == value2 else ''


@register.simple_tag
def type_label(data):
    if data.type:
        label = lambda x: f"<span class='ui label' > {x} </span>"
        types = [label(t) for t in data.type.split(',')]
        return mark_safe(''.join(types))
    return ""


@register.simple_tag
def img(obj):
    """
    Returns the image associated with the object or a placeholder
    """
    if obj.image:
        return obj.image.url
    else:
        return urllib.parse.urljoin(settings.STATIC_URL, "images/placeholder.png")


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.inclusion_tag('widgets/recipe_form.html')
def recipe_form(form):
    """
    Renders a recipe form.
    """
    return dict(form=form)


@register.inclusion_tag('parts/recipe_details.html', takes_context=True)
def recipe_details(context, recipe, include_copy=True):
    user = context['request'].user

    return dict(user=user, recipe=recipe, project=recipe.project, include_copy=include_copy)


@register.filter
def writable(project, user):
    """
    Check if user has write access to the project.
    """
    return auth.is_writable(user=user, project=project)


@register.simple_tag
def image_field(default=''):
    if default:
        image_field = forms.ImageField(required=False, default=default)
    else:
        image_field = forms.ImageField(required=False)
    image_field.widget.attrs.update({'id': 'image'})
    placeholder = urllib.parse.urljoin(settings.STATIC_URL, "images/placeholder.png")
    image_widget = image_field.widget.render('image', value=placeholder)

    return mark_safe(image_widget)


@register.inclusion_tag('widgets/created_by.html')
def created_by(date, user=None, prefix="updated"):
    """
    Renders a created by link
    """
    return dict(date=date, user=user, prefix=prefix)

@register.inclusion_tag('widgets/recipe_clone_message.html')
def recipe_clone_message(recipe):
    """
    Renders the recipe clone message.
    """
    return dict(recipe=recipe)


@register.inclusion_tag('widgets/loading_img.html')
def job_img(job):
    return dict(job=job)


@register.inclusion_tag('widgets/access_form.html')
def access_form(project, user, extra_class=''):
    """
    Generates an access form.
    """

    return dict(project=project, user=user, extra_class=extra_class)


@register.filter
def get_access_label(user, project):
    access = Access.objects.filter(user=user, project=project).select_related('project').first()

    access = access or Access(access=Access.NO_ACCESS, user=user, project=project)

    label = access.get_access_display()
    return label


@register.filter
def get_access(user, project):
    access = Access.objects.filter(user=user, project=project).first()
    access = access or Access(access=Access.NO_ACCESS, user=user, project=project)
    return access


@register.inclusion_tag('widgets/job_elapsed.html')
def job_minutes(job, view=False):
    check_back = ''
    # Add a tag to check a state change every ~5 seconds and update tag
    if job.state in [Job.SPOOLED, Job.RUNNING, Job.QUEUED]:
        check_back = 'check_back'

    return dict(job=job, check_back=check_back, view=view)


@register.simple_tag
def size_label(data):
    """
    Returns a label for data sizes.
    """

    size = f"{defaultfilters.filesizeformat(data.size)}"
    return mark_safe(f"<span class='ui mini label'>{size}</span>")


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


def now():
    return datetime.utcnow().replace(tzinfo=utc)

def pluralize(value, word):
    if value > 1:
        return "%d %ss" % (value, word)
    else:
        return "%d %s" % (value, word)

@register.filter
def time_ago(date):
    # Rare bug. TODO: Need to investigate why this can happen.
    if not date:
        return ''
    delta = now() - date
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
