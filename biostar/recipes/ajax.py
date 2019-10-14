from functools import wraps, partial
import logging
import hjson
from ratelimit.decorators import ratelimit

from django.template import loader
from django.http import JsonResponse
from django.utils.decorators import available_attrs
from django.template import loader
from .const import *
from biostar.recipes.models import Job, Analysis
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
                #1/0
                return ajax_error(f'{self.method} method must be used.')
            #1/0
            if not request.user.is_authenticated:
                #1/0
                return ajax_error('You must be logged in.')
            #1/0
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


@ratelimit(key='ip', rate='50/h')
@ratelimit(key='ip', rate='10/m')
@ajax_error_wrapper(method="GET")
def preview_json(request):

    # Get the recipe
    recipe_uid = request.GET.get('uid')

    recipe = Analysis.objects.filter(uid=recipe_uid).first()

    if not recipe:
        return ajax_error(msg="Recipe does not exist.")

    json_text = request.GET.get('json_text', recipe.json_text)
    json_data = hjson.loads(json_text)
    # Render the recipe interface
    interface = RecipeInterface(request=request, analysis=recipe, json_data=json_data,
                                initial=dict(name=recipe.name), add_captcha=False)

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
        opts.update(dict(source=source))
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
        dropdown=(DROPDOWN, None, [('Choices 1', 'Choices 1'), ('Choices 2', 'Choices 2')], 'Choices 1',
                  'Dropdown Field Label'),
        integer=(INTEGER, None, [], 100, 'Integer Field Label'),
        textbox=(TEXTBOX, None, [], "Sample text", 'Text box Field Label'),
        float=(FLOAT, None, [], 0.5, 'Float Field Label'),
        checkbox=(CHECKBOX, None, [], True, 'Checkbox Field Label'),
    )
    display, source, choices, value, label = display_map.get(display_type, ('textbox', []))
    #print(display, choices)
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
