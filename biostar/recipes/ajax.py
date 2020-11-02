from functools import wraps, partial
import logging

from django.contrib import messages
import toml
from ratelimit.decorators import ratelimit
from django.views.decorators.csrf import ensure_csrf_cookie
from django.http import HttpResponse
from django.template import Template, Context
from django.http import JsonResponse
from django.template import loader
from django.conf import settings
from biostar.accounts.models import User
from biostar.recipes.const import *
from biostar.recipes.models import Job, Analysis, Data, Project, MAX_TEXT_LEN, Access
from biostar.recipes.forms import RecipeInterface
from biostar.recipes import auth, util
from django.shortcuts import render, redirect, reverse

from biostar.recipes.forms import RecipeForm

logger = logging.getLogger("engine")

JOB_COLORS = {Job.SPOOLED: "spooled",
              Job.ERROR: "errored", Job.QUEUED: "queued",
              Job.RUNNING: "running", Job.COMPLETED: "completed"
              }

RATELIMIT_KEY = settings.RATELIMIT_KEY

def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return JsonResponse(payload)


ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')

MIN_TITLE_CHARS = 10
MAX_TITLE_CHARS = 180
MAX_SNIPPET = 180
MAX_HELP = 180
MAX_SNIPPETS_CATEGORIES = 10
MAX_SNIPPETS_PER_CATEGORY = 30

MAX_TAGS = 5


class ajax_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, method, login_required=True):
        self.method = method
        self.login_required = login_required

    def __call__(self, func, *args, **kwargs):

        @wraps(func)
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')
            if request.user.is_anonymous and self.login_required:
                return ajax_error('You must be logged in.')
            return func(request, *args, **kwargs)

        return _ajax_view


@ajax_error_wrapper(method="POST", login_required=True)
def ajax_edit(request, id):
    """
    Edit recipes. Uses recipe primary key as access!
    """

    # The user performing the request.
    user = request.user

    # The recipe that needs to be edited.
    recipe = Analysis.objects.filter(id=id).first()

    if not auth.writeable_recipe(source=recipe, user=user):
        message = str("Recipe is not writable by current user")
        return ajax_error(msg=message)

    # The project of for the recipe.
    project = recipe.project
    form = RecipeForm(data=request.POST, instance=recipe, files=request.FILES, user=user,
                      project=project)

    if form.is_valid():
        form.save()
        return ajax_success(msg=f"Recipe saved.")
    else:
        message = str(form.errors)
        return ajax_error(msg=message)


def check_job(request, uid):
    job = Job.objects.filter(uid=uid).first()

    check_back = 'check_back' if job.state in [Job.SPOOLED, Job.RUNNING] else ''
    stdout_path = os.path.join(job.path, settings.JOB_STDOUT)
    stderr_path = os.path.join(job.path, settings.JOB_STDERR)

    try:
        state_changed = int(request.GET.get('state')) != job.state
    except Exception as exc:
        logger.error(f'Error checking job:{exc}')
        state_changed = False

    if os.path.exists(stdout_path) and os.path.exists(stderr_path):
        stdout = open(stdout_path, 'r').read()
        stderr = open(stderr_path, 'r').read()
        Job.objects.filter(uid=job.uid).update(stderr_log=stderr, stdout_log=stdout)
    else:
        stdout = stderr = None

    # Render the updated image icon
    redir = job.url() if job.is_finished() else ""
    context = dict(job=job)
    tmpl = loader.get_template('widgets/loading_img.html')
    image_tmpl = tmpl.render(context=context)

    # Render the updated elapsed runtime, status, etc..
    context = dict(job=job, check_back=check_back)
    tmpl = loader.get_template('widgets/job_elapsed.html')
    template = tmpl.render(context=context)

    return ajax_success(msg='success', redir=redir,
                        html=template, state=job.get_state_display(), is_running=job.is_running(),
                        stdout=stdout, stderr=stderr, job_color=auth.job_color(job),
                        state_changed=state_changed, img_tmpl=image_tmpl)


def check_size(fobj, maxsize=0.3):
    # maxsize in megabytes!

    try:
        if fobj and fobj.size > maxsize * 1024 * 1024.0:
            curr_size = fobj.size / 1024 / 1024.0
            msg = f"File too large: {curr_size:0.1f}MB should be < {maxsize:0.1f}MB"
            return msg, False
    except Exception as exc:
        return f"File size validation error: {exc}", False

    return "Valid", True


@ajax_error_wrapper(method="POST", login_required=False)
def preview_json(request):
    # Get the recipe
    uid = request.POST.get('recipe')
    text = request.POST.get('toml', '')
    recipe = Analysis.objects.filter(id=uid).first()

    if not recipe:
        return ajax_error(msg=f"Project does not exist")

    try:
        json_data = toml.loads(text)
    except toml.decoder.TomlDecodeError as exc:
        msg = util.toml_error(exp_msg=exc, text=text)
        return ajax_error(msg=f"{msg}")

    project = recipe.project
    # Render the recipe interface
    interface = RecipeInterface(request=request, json_data=json_data, project=project,
                                add_captcha=False)

    tmpl = loader.get_template('widgets/recipe_form.html')
    context = dict(form=interface, focus=True)
    template = tmpl.render(context=context)

    return ajax_success(msg="Recipe json", html=template)


def drop(request, klass=None):
    """
    Order objects in a list using drag and drop. 'rank' attribute is required .
    """
    user = request.user
    # Current object uid being moved/dropped
    source = int(request.POST.get("source_id", 0))
    source = klass.objects.filter(pk=source).first()

    # The object we intend to move it under
    top = int(request.POST.get("parent_id", 0))
    top = klass.objects.filter(pk=top).first()

    # Next object after the 'source'
    bottom = int(request.POST.get("next_id", 0))
    bottom = klass.objects.filter(pk=bottom).first()

    project = source.project
    # Check if the user has write access to source before moving.
    if not auth.is_writable(user=user, project=project):
        ajax_error(msg="You need write access to move objects.")

    # Compute and update the source with a new rank.
    maxrank = klass.objects.order_by('-rank').first().rank
    source.rank = auth.compute_rank(source=source, top=top, bottom=bottom, klass=klass, maxrank=maxrank)
    klass.objects.filter(pk=source.pk).update(rank=source.rank)

    return ajax_success(msg="Successfully moved")


@ajax_error_wrapper(method="POST", login_required=True)
def drop_project(request):
    return drop(request=request, klass=Project)


@ajax_error_wrapper(method="POST", login_required=True)
def drop_recipe(request):
    return drop(request=request, klass=Analysis)


@ajax_error_wrapper(method="POST", login_required=True)
def toggle_delete(request):
    """
    Delete an object.
    """

    type_map = dict(job=Job, data=Data, recipe=Analysis)
    uid = request.POST.get('uid', "")
    obj_type = request.POST.get('type', '')

    obj_model = type_map.get(obj_type)

    if not obj_model:
        return ajax_error(f"Invalid data type:{obj_type}")

    obj = obj_model.objects.filter(uid=uid).first()

    if not obj:
        return ajax_error("Object does not exists.")

    access = auth.is_writable(user=request.user, project=obj.project)

    # Toggle the delete state if the user has write access
    if access:
        auth.delete_object(obj=obj, request=request)
        # Re-set project counts
        obj.project.set_counts(save=True)
        counts = obj_model.objects.filter(project=obj.project, deleted=False).count()
        return ajax_success(msg='Toggled delete', counts=counts)

    return ajax_error("Invalid action")


@ratelimit(key=RATELIMIT_KEY, rate='10/m')
@ajax_error_wrapper(method="POST", login_required=True)
def manage_access(request):
    access_map = dict(none=Access.NO_ACCESS, read=Access.READ_ACCESS,
                      write=Access.WRITE_ACCESS, share=Access.SHARE_ACCESS)

    # Get the current user, project and access
    user_id = request.POST.get('user_id', '')
    project_uid = request.POST.get('project_uid', '')
    access_str = request.POST.get('access', '')

    user = User.objects.filter(id=user_id).first()
    new_access = access_map.get(access_str)
    project = Project.objects.filter(uid=project_uid).first()
    is_writable = auth.is_writable(user=request.user, project=project)

    # Validate submitted values.
    if not project:
        return ajax_error("Project does not exist.")
    if not user:
        return ajax_error("User does not exist.")
    if not new_access:
        return ajax_error(f"Invalid access option: {access_str}")
    if user == request.user:
        return ajax_error("Can not change your own access")
    if user == project.owner:
        return ajax_error("Can not change the project owner's access")
    if not is_writable:
        return ajax_error("You need write access to manage access.")

    # Check current user access.
    access = Access.objects.filter(user=user, project=project).first()

    # Update existing access
    if access:
        Access.objects.filter(id=access.id).update(access=new_access)
    # Create a new access object
    else:
        Access.objects.create(user=user, project=project, access=new_access)

    no_access = new_access == Access.NO_ACCESS

    return ajax_success("Changed access.", no_access=no_access)

@ensure_csrf_cookie
@ajax_error_wrapper(method="POST", login_required=True)
def copy_file(request):
    """
    Add file from import root directory into clipboard.
    """
    # Assumes incoming path is a full path
    path = request.POST.get('path')
    user = request.user

    # Project uid to check access
    uid = request.POST.get('uid') or 0
    project = Project.objects.filter(id=uid).first()

    if uid and not project:
        return ajax_error(msg="Project does not exist.")

    if project and not auth.is_readable(user=user, obj=project):
        return ajax_error(msg="You do not have access to copy this file")

    if not os.path.exists(path):
        return ajax_error(msg="File path does not exist.")

    if path.startswith(settings.IMPORT_ROOT_DIR) and not user.profile.trusted:
        return ajax_error(msg="Only trusted users can copy from this directory.")

    copied = auth.copy_file(request=request, fullpath=path)
    return ajax_success(msg=f"{len(copied)} files copied.")


@ensure_csrf_cookie
@ajax_error_wrapper(method="POST", login_required=True)
def copy_object(request):
    """
    Add object uid or file path to sessions clipboard.
    """
    # Map the query parameter to a clipboard and database model.
    mapper = {"data": (Data, COPIED_DATA),  "job": (Job, COPIED_RESULTS), "recipe": (Analysis, COPIED_RECIPES)}
    uid = request.POST.get('uid', '')
    oid = request.POST.get('id', 0)

    # Get the clipboard to copy to
    clipboard = request.POST.get(settings.CLIPBOARD_NAME)

    klass, board = mapper.get(clipboard, (None, None))
    if uid:
        obj = klass.objects.filter(uid=uid).first()
    else:
        obj = klass.objects.filter(id=oid).first()

    if not (klass and board and obj):
        return ajax_error("Object or board does not exist.")

    is_readable = auth.is_readable(user=request.user, obj=obj)

    if not is_readable:
        return ajax_error('You do not have access to copy this object.')

    # Return current clipboard contents
    copied = auth.copy_uid(request=request, uid=obj.uid, board=board)

    return ajax_success(f"Copied. Clipboard contains :{len(copied)} objects.")


@ajax_error_wrapper(method="POST", login_required=True)
def ajax_clear_clipboard(request):

    # Clear the clipboard
    auth.clear(request=request)

    return ajax_success(msg="Cleared clipboard. ", html="")


@ajax_error_wrapper(method="POST", login_required=True)
def ajax_move(request):

    pid = request.POST.get("id", 0)
    user = request.user
    project = Project.objects.filter(id=pid).first()

    # Get the board.
    board = auth.recent_clipboard(request=request)
    key, vals = board
    next_url = auth.resolve_paste_url(key=key, project=project)

    count = len(vals)

    if not project:
        return ajax_error(msg="Project does not exist.")

    if not auth.is_writable(user=user, project=project):
        return ajax_error(msg="You do not have access to move here.")

    # Move objects in clipboard to given project.
    auth.move(uids=vals, project=project, otype=key, user=user)

    # Clear the clipboard after moving.
    auth.clear(request=request)

    return ajax_success(msg=f"Moved {count} items into project.", redirect=next_url)


@ajax_error_wrapper(method="POST", login_required=True)
def ajax_paste(request):
    """
    Paste the most recent
    """
    pid = request.POST.get("id", 0)
    user = request.user
    project = Project.objects.filter(id=pid).first()

    # Get the board.
    board = auth.recent_clipboard(request=request)
    key, vals = board
    count = len(vals)

    if not project:
        return ajax_error(msg="Project does not exist.")

    if not auth.is_writable(user=user, project=project):
        return ajax_error(msg="You do not have access to paste here.")

    if not count:
        return ajax_error(msg="Clipboard is empty")

    # The target of this action is to clone.
    clone = request.POST.get('target') == CLONED_RECIPES

    # Paste the clipboard item into the project
    auth.paste(board=board, user=user, project=project, clone=clone)

    # Resolve the redirect url.
    next_url = auth.resolve_paste_url(key=key, project=project)

    # Clear the clipboard after pasting.
    auth.clear(request=request)

    return ajax_success(msg=f"Pasted {count} items into project.", redirect=next_url)


@ajax_error_wrapper(method="GET", login_required=False)
def render_plugins(request):

    # Plugin to apply.
    plugin = request.GET.get('plugin')
    # Filename to apply plugin to.
    fname = request.GET.get('fname')
    html = ''
    if plugin == 'jupyter':
        html = open(fname, 'r').read()

    return ajax_success(msg="Processed successfully", html=html)


@ensure_csrf_cookie
@ajax_error_wrapper(method="POST", login_required=True)
def ajax_clipboard(request):
    """
    Displays the most recent object in clipboard.
    """

    # Get the most current item in the clipboard.
    pid = request.POST.get("id", 0)
    user = request.user

    project = Project.objects.filter(id=pid).first()
    board = auth.recent_clipboard(request=request)
    key, vals = board
    count = len(vals)

    if project and auth.is_readable(user=user, obj=project) and count:
        # Load items into clipboard
        tmpl = loader.get_template('widgets/clipboard.html')
        movable = key in [COPIED_RECIPES, COPIED_DATA]
        context = dict(count=count, board=key, is_recipe=key == COPIED_RECIPES, movable=movable)
        template = tmpl.render(context=context)
    else:
        template = ''

    return ajax_success(html=template, msg="Refreshed clipboard")


def field_render(request):
    """
    Renders and returns an HTML field based on the incoming toml data.
    """

    demo = """
    [demo]
    label = "Demo"
    display = "INTEGER"
    value = 10000000
    """
    text = request.POST.get("toml", demo)

    # Get the project uid
    uid = request.POST.get('recipe')
    recipe = Analysis.objects.filter(id=uid).first()
    try:
        data = toml.loads(text)
    except Exception as exc:
        return HttpResponse(exc)

    form = RecipeInterface(request, json_data=data, analysis=recipe)

    context = dict(form=form)

    return render(request, "parts/form_field_render.html", context=context)

