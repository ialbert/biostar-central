from functools import wraps, partial
import logging
import hjson, json
from ratelimit.decorators import ratelimit

from django.shortcuts import reverse
from django.template import loader
from django.template import Template, Context
from django.http import JsonResponse
from django.utils.decorators import available_attrs
from django.template import loader
from django.conf import settings
from biostar.recipes.const import *
from biostar.recipes.models import Job, Analysis, Snippet, SnippetType, Project, MAX_TEXT_LEN
from biostar.recipes.forms import RecipeInterface
from biostar.recipes import auth

logger = logging.getLogger("engine")

JOB_COLORS = {Job.SPOOLED: "spooled",
              Job.ERROR: "errored", Job.QUEUED: "queued",
              Job.RUNNING: "running", Job.COMPLETED: "completed"
              }

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

        @wraps(func, assigned=available_attrs(func))
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')
            if request.user.is_anonymous and self.login_required:
                return ajax_error('You must be logged in.')
            return func(request, *args, **kwargs)

        return _ajax_view


def check_job(request, uid):

    job = Job.objects.filter(uid=uid).first()
    check_back = 'check_back' if job.state in [Job.SPOOLED, Job.RUNNING] else ''

    try:
        state_changed = int(request.GET.get('state')) != job.state
    except Exception as exc:
        logger.error(f'Error checking job:{exc}')
        state_changed = False

    context = dict(job=job, check_back=check_back)
    tmpl = loader.get_template('widgets/job_elapsed.html')
    template = tmpl.render(context=context)

    return ajax_success(msg='success', html=template, state=job.get_state_display(), state_changed=state_changed)


@ajax_error_wrapper(method="POST", login_required=False)
def snippet_code(request):

    command_uid = request.POST.get('command', '')
    current_code = request.POST.get('template', '')

    cmd = Snippet.objects.filter(uid=command_uid).first()

    if not cmd:
        return ajax_error(msg='Command not found.')

    command = cmd.command
    comment = f'# { cmd.help_text }' if cmd.help_text else ' '
    code = current_code + '\n' + comment + '\n' + command

    tmpl = loader.get_template('widgets/template_field.html')
    context = dict(template=code, scroll_to_bottom=True)
    template_field = tmpl.render(context=context)

    return ajax_success(code=code, msg="Rendered the template", html=template_field)


#@ratelimit(key='ip', rate='50/h')
#@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST", login_required=False)
def snippet_form(request):

    is_category = request.POST.get("is_category", 0)
    is_category = bool(int(is_category))

    type_name = request.POST.get('type_name')
    type_uid = request.POST.get('type_uid')
    snippet_uid = request.POST.get('snippet_uid', '')
    snippet = request.POST.get('snippet', '')
    help_text = request.POST.get('help_text', '')

    tmpl = loader.get_template('widgets/snippet_form.html')

    context = dict(is_category=is_category, type_name=type_name, snippet=snippet, help_text=help_text,
                   type_uid=type_uid, snippet_uid=snippet_uid)
    cmd_form = tmpl.render(context=context)

    return ajax_success(html=cmd_form, msg="Rendered form")


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


#@ratelimit(key='ip', rate='50/h')
#@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def create_snippet_type(request):
    name = request.POST.get('name', '')
    image = request.FILES.get('image', '')

    if not name:
        return ajax_error(msg="Name is required.")

    user_types = SnippetType.objects.filter(owner=request.user)

    if user_types.count() >= MAX_SNIPPETS_CATEGORIES and not request.user.is_superuser:
        return ajax_error(msg="Maximum amount of snippets reached.")

    # Get the type of code this is: Bash, r, Matlab, etc...
    cmd_type = SnippetType.objects.create(name=name, owner=request.user)

    if image:
        msg, valid = check_size(fobj=image)
        if not valid:
            return ajax_error(msg=msg)
        cmd_type.image = image
        cmd_type.save()

    tmpl = loader.get_template('widgets/snippet_type.html')
    context = dict(type=cmd_type, request=request)
    new_type = tmpl.render(context=context)

    return ajax_success(msg="Created snippet", html=new_type)


#@ratelimit(key='ip', rate='50/h')
#@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def create_snippet(request):

    snippet = request.POST.get('snippet', '')
    help_text = request.POST.get('help_text', '')
    type_uid = request.POST.get('type_uid')
    snippet_uid = request.POST.get('snippet_uid', '')

    if not (help_text and snippet):
        return ajax_error(msg="Snippet and help text are required.")

    # Get the type of code this this: Bash, r, Matlab, etc...
    cmd_type = SnippetType.objects.filter(uid=type_uid).first()
    if not cmd_type:
        return ajax_error(msg=f"Snippet does not have a valid type:{cmd_type}")

    snippets = Snippet.objects.filter(type__owner=request.user, type=cmd_type)

    if snippets.count() >= MAX_SNIPPETS_PER_CATEGORY:
        return ajax_error(msg="Maximum number of snippets reached.")

    if len(snippet) >= MAX_SNIPPET or len(help_text) >= MAX_HELP:
        msg = "Snippet" if len(snippet) >= MAX_TEXT_LEN else "Help Text"
        return ajax_error(msg=msg + " input is too long.")

    # Get existing snippet for edit
    if snippet_uid:
        snippet_obj = Snippet.objects.filter(uid=snippet_uid).first()
        if not snippet:
            return ajax_error(msg=f"Editing error: snippet id does not exist {snippet_uid}.")
        # Update the snippet and help_text
        Snippet.objects.filter(uid=snippet_obj.uid).update(help_text=help_text, command=snippet)
        # Re-fetch the snippet after the edit.
        snippet = Snippet.objects.filter(uid=snippet_obj.uid).first()
    else:
        snippet = Snippet.objects.create(command=snippet, type=cmd_type, help_text=help_text, owner=request.user)

    # Load snippet into template
    tmpl = loader.get_template('widgets/snippet.html')
    context = dict(snippet=snippet)
    created_form = tmpl.render(context=context)

    return ajax_success(msg="Created snippet", html=created_form)


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def delete_snippet(request):

    snippet_uid = request.POST.get('snippet_uid', '')
    # Get the snippet
    snippet = Snippet.objects.filter(uid=snippet_uid).first()

    if request.user != snippet.owner:
        return ajax_error(msg="Only owners or superusers can delete their code snippets.")

    return


@ajax_error_wrapper(method="POST", login_required=False)
def preview_template(request):

    source_template = request.POST.get('template', '# Code goes here')
    source_json = request.POST.get('json_text', '{}')
    name = request.POST.get('name', 'Name')
    project_uid = request.POST.get('project_uid')
    project = Project.objects.filter(uid=project_uid).first()

    try:
        source_json = hjson.loads(source_json)
        # Fill json information by name for the preview.
        source_json = auth.fill_data_by_name(project=project, json_data=source_json)
        # Fill in the script with json data.
        context = Context(source_json)
        script_template = Template(source_template)
        script = script_template.render(context)
    except Exception as exc:
        return ajax_error(f"Error rendering code: {exc}")

    # Load the html containing the script
    tmpl = loader.get_template('widgets/preview_template.html')
    context = dict(script=script, name=name)
    template = tmpl.render(context=context)

    return ajax_success(html=template, msg="Rendered script")


@ajax_error_wrapper(method="POST", login_required=False)
def preview_json(request):

    # Get the recipe
    recipe_name = request.POST.get('name')
    project_uid = request.POST.get('project_uid')

    json_text = request.POST.get('json_text', '{}')
    try:
        json_data = hjson.loads(json_text)
    except Exception as exc:
        return ajax_error(msg=f"{exc}")

    project = Project.objects.filter(uid=project_uid).first()
    if not project:
        return ajax_error(msg=f"Project does not exist")

    # Render the recipe interface
    interface = RecipeInterface(request=request, json_data=json_data, project=project,
                                initial=dict(name=recipe_name), add_captcha=False)

    tmpl = loader.get_template('widgets/recipe_form.html')
    context = dict(form=interface, focus=True)
    template = tmpl.render(context=context)

    return ajax_success(msg="Recipe json", html=template)


def get_display_dict(display_type):
    mapping = dict(radio=RADIO, integer=INTEGER, textbox=TEXTBOX,
                   float=FLOAT, checkbox=CHECKBOX, dropdown=DROPDOWN,
                   upload=UPLOAD)

    display = mapping.get(display_type)

    if display_type == 'data':
        return dict(label='Data Field Label',
                    source='PROJECT',
                    type='DATA',
                    help='Pick data from this project to analyze.')
    if display == RADIO:
        return dict(label='Radio Field Label',
                    display=RADIO, help='Choose an option.',
                    choices=[(1, 'Option 1'), (2, 'Option 2')], value=2)
    if display == INTEGER:
        return dict(label='Integer Field Label',
                    display=INTEGER,
                    help='Enter an integer between -100 and 100.',
                    range=[-100, 100], value=0)
    if display == TEXTBOX:
        return dict(label='Text box Field Label', display=TEXTBOX,
                    help='Enter text.',
                    value='text')
    if display == FLOAT:
        return dict(label='Float Field Label',
                    help='Enter a float, decimal number, between -100.0 and 100.0.',
                    display=FLOAT, range=[-100.0, 100.0],
                    value=0.5)
    if display == CHECKBOX:
        return dict(label='Checkbox Field Label',
                    help="Check the box for 'yes'. ",
                    display=CHECKBOX, value=True)
    if display == DROPDOWN:
        return dict(label='Dropdown Field Label',
                    display=DROPDOWN,
                    help="Pick an option from a dropdown.",
                    choices=[('1', 'Choices 1'), ('2', 'Choices 2')],
                    value='1')
    if display == UPLOAD:
        return dict(label='Upload a file',
                    display=UPLOAD,
                    help="Upload a file to analyze")

    return dict()


@ajax_error_wrapper(method="POST", login_required=False)
def add_to_interface(request):
    # Returns a recipe interface field json

    display_type = request.POST.get('display_types', '')
    json_text = request.POST.get('json_text', '')

    display_dict = get_display_dict(display_type=display_type)

    json_data = hjson.loads(json_text)
    field_name = display_type
    count = 0
    # Check if the field name exists
    while field_name in json_data:
        field_name = display_type + f'{count}'
        count += 1

    new_field = {field_name: display_dict}
    json_data.update(new_field)
    new_json = hjson.dumps(json_data)

    tmpl = loader.get_template('widgets/json_field.html')
    context = dict(json_text=new_json, focus=True)
    json_field = tmpl.render(context=context)

    return ajax_success(html=json_field, json_text=new_json, msg="Rendered json")


@ajax_error_wrapper(method="POST")
def file_copy(request):
    """
    Add file into clipboard.
    """
    path = request.POST.get('path')
    fullpath = os.path.abspath(os.path.join(settings.IMPORT_ROOT_DIR, path))

    if not os.path.exists(fullpath):
        return ajax_error(msg="File path does not exist.")

    copied = auth.copy_file(request=request, fullpath=fullpath)

    return ajax_success(msg=f"{len(copied)} files copied.")


@ajax_error_wrapper(method="POST", login_required=True)
def toggle_delete(request):
    """
    Delete a job
    """
    job_uid = request.POST.get('job_uid', "")
    job = Job.objects.filter(uid=job_uid).first()

    if not job:
        return ajax_error("Job does not exists.")

    access = auth.is_writable(user=request.user, project=job.project)

    # Toggle the delete state if the user has write access
    if access:
        deleted = auth.delete_object(obj=job, request=request)
        msg_prefix = "Deleted" if deleted else "Restored"
        counts = Job.objects.filter(project=job.project, deleted=False).count()
        return ajax_success(f"{msg_prefix} {job.name}", counts=counts)

    return ajax_error("Invalid action")


@ajax_error_wrapper(method="POST", login_required=True)
def copy_object(request):
    """
    Add object uid or file path to sessions clibboard.
    """

    object_uid = request.POST.get('uid', '')
    project_uid = request.POST.get('project_uid', '')
    clipboard = request.POST.get('clipboard')

    project = Project.objects.filter(uid=project_uid).first()
    if not project:

        return ajax_error("Project does not exist.")

    is_readable = auth.is_readable(user=request.user, project=project)

    if not is_readable:
        return ajax_error('You do not have access to copy this object.')

    # Return current clipboard contents
    copied_uids = auth.copy_uid(request=request, uid=object_uid, board=clipboard)
    clipboard = request.session.get(settings.CLIPBOARD_NAME, {})

    print(clipboard)

    return ajax_success(f"Copied. Clipboard contains :{len(copied_uids)} objects.")


def add_variables(request):

    # Get the most recent template and json.
    json_text = request.POST.get('json_text', '')
    template = request.POST.get('template', '')

    json_data = hjson.loads(json_text)

    # Create a set with all template variables
    all_vars = {"{{ " + f"{v}.value" + "}}" for v in json_data.keys()}

    all_vars = '\n'.join(all_vars)

    # Insert variables at the beginning of the template
    new_template = template + '\n' + all_vars

    tmpl = loader.get_template('widgets/template_field.html')
    context = dict(template=new_template, scroll_to_bottom=False, focus=True)
    template_field = tmpl.render(context=context)

    return ajax_success(msg="Added variables to template", html=template_field, code=new_template)