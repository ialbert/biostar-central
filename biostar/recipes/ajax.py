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
from .const import *
from biostar.recipes.models import Job, Analysis, Snippet, SnippetType,Project, MAX_TEXT_LEN
from biostar.recipes.forms import RecipeInterface

logger = logging.getLogger("engine")


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

    def __init__(self, method):
        self.method = method

    def __call__(self, func, *args, **kwargs):

        @wraps(func, assigned=available_attrs(func))
        def _ajax_view(request, *args, **kwargs):

            if request.method != self.method:
                return ajax_error(f'{self.method} method must be used.')
            if not request.user.is_authenticated:
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


#@ratelimit(key='ip', rate='50/h')
#@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
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
    context = dict(template=code)
    template_field = tmpl.render(context=context)

    return ajax_success(code=code, msg="Rendered the template", html=template_field)


#@ratelimit(key='ip', rate='50/h')
#@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def snippet_form(request):

    is_top = request.POST.get("is_top", "false")

    is_top = True if is_top == 'true' else False
    type_name = request.POST.get('type_name')
    type_uid = request.POST.get('type_uid')
    snippet_uid = request.POST.get('snippet_uid', '')
    snippet = request.POST.get('snippet', '')
    help_text = request.POST.get('help_text', '')

    tmpl = loader.get_template('widgets/snippet_form.html')

    context = dict(is_top=is_top, type_name=type_name, snippet=snippet, help_text=help_text,
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


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
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


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
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

    cmd = snippets.filter(command=snippet).first()
    # Avoid duplicates commands being created
    if cmd and not snippet_uid:
        return ajax_error(msg=f"Similar snippet for {cmd_type.name} already exists: {snippet}.")

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


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def preview_template(request):

    source_template = request.POST.get('template', '# Code goes here')
    source_json = request.POST.get('json_text', '{}')
    name = request.POST.get('name', 'Name')
    recipe_uid = request.POST.get('uid')
    source_json = hjson.loads(source_json)

    # Load the recipe JSON into the template
    context = Context(source_json)
    script_template = Template(source_template)
    script = script_template.render(context)

    recipe = Analysis.objects.filter(uid=recipe_uid).first()
    download_url = reverse('recipe_download', kwargs=dict(uid=recipe_uid)) if recipe else '#template_field'

    # Load the html containing the script
    tmpl = loader.get_template('widgets/preview_download.html')
    context = dict(script=script, name=name, download_url=download_url)
    template = tmpl.render(context=context)

    return ajax_success(html=template, msg="Rendered script")


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
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
    context = dict(form=interface)
    template = tmpl.render(context=context)

    return ajax_success(msg="Recipe json", html=template)


def field_specs(display_type='', source=None, choices='', value='', label="label"):

    opts = dict(label=label or f'{display_type.lower()}')
    if display_type:
        opts.update(dict(display=display_type))
    if choices:
        opts.update(dict(choices=choices))
    if source:
        opts.update(dict(source=source, type="DATA"))
    if value:
        opts.update(dict(value=value))

    return opts


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="POST")
def recipe_field(request):
    # Returns a recipe interface field json

    display_type = request.POST.get('display_types', '')
    json_text = request.POST.get('json_text', '')

    display_map = dict(
        radio=(RADIO, None, [(1, 'Option 1'), (2, 'Option 2')], 2, 'Radio Field Label'),
        data=('', 'PROJECT', [], '', 'Data Field Label'),
        integer=(INTEGER, None, [], 100, 'Integer Field Label'),
        textbox=(TEXTBOX, None, [], "Sample text", 'Text box Field Label'),
        float=(FLOAT, None, [], 0.5, 'Float Field Label'),
        checkbox=(CHECKBOX, None, [], True, 'Checkbox Field Label'),
        dropdown=(DROPDOWN, None, [('Choices 1', 'Choices 1'), ('Choices 2', 'Choices 2')], 'Choices 1',
                  'Dropdown Field Label'),
    )
    display, source, choices, value, label = display_map.get(display_type, ('textbox', []))
    # Return the field specs
    specs = field_specs(display_type=display, source=source, choices=choices, value=value, label=label)

    json_data = hjson.loads(json_text)
    field_name = display_type
    count = 0
    # Check if the field name exists
    while field_name in json_data:
        field_name = display_type + f'{count}'
        count += 1

    new_field = {field_name: specs}
    new_field.update(json_data)
    new_json = hjson.dumps(new_field)

    tmpl = loader.get_template('widgets/json_field.html')
    context = dict(json_text=new_json)
    json_field = tmpl.render(context=context)

    return ajax_success(html=json_field, json_text=new_json, msg="Rendered json")


def add_variables(request):

    json_text = request.POST.get('json_text', '')
    template = request.POST.get('template', '')

    vars = hjson.loads(json_text).keys()
    vars = ["{{ " + f"{v}.value" + "}}" for v in vars]

    vars = '\n'.join(vars)

    new_template = vars + "\n" + template

    tmpl = loader.get_template('widgets/template_field.html')
    context = dict(template=new_template)
    template_field = tmpl.render(context=context)

    return ajax_success(msg="Added variables to template", html=template_field, code=new_template)