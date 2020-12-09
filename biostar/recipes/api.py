import toml
import json
import logging
import os, urllib
import base64
from functools import wraps
from django.conf import settings
from django.http import HttpResponse
from ratelimit.decorators import ratelimit
from django.views.decorators.csrf import csrf_exempt
from biostar.accounts.models import User
from biostar.recipes.models import Analysis, Project, Data, image_path, Access
from biostar.recipes import util, auth
from biostar.recipes.decorators import token_access, require_api_key

logger = logging.getLogger("engine")

RATELIMIT_KEY = settings.RATELIMIT_KEY

# Maximum file size to be sent and received via api.
MAX_FILE_SIZE = 100


class api_error_wrapper:
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, methods=['GET']):
        self.methods = methods

    def __call__(self, func, *args, **kwargs):
        @wraps(func)
        def _ajax_view(request, *args, **kwargs):
            if request.method not in self.methods:
                return HttpResponse(content=f'{self.methods} method must be used.')

            return func(request, *args, **kwargs)

        return _ajax_view


def get_thumbnail():
    return os.path.join(settings.STATIC_ROOT, "images", "placeholder.png")


def change_image(obj, file_object=None):
    if not obj:
        return get_thumbnail()

    obj.image.save(name=get_thumbnail(), content=file_object)

    return obj.image.path


def encode_image(img):
    """
    Base 64 encoding of an image field that may be missing.
    """

    # The image is not filled in.
    if not img:
        return ''

    # The path is incorrect.
    if not os.path.isfile(img.path):
        return ''

    # The binary representation of the data.
    data = open(img.path, 'rb').read()

    # Convert image to base64 ASCII string
    text = base64.b64encode(data).decode("ascii")

    return text


def encode_project(project, show_image=False):
    recipes = dict()
    store = dict(
        uid=project.uid,
        name=project.name,
        text=project.text,
        date=str(project.date),
        privacy=project.privacy,
        image=encode_image(project.image) if show_image else '',
        recipes=recipes,
    )
    return store


def encode_recipe(recipe, show_image=False):
    store = dict(
        uid=recipe.uid,
        name=recipe.name,
        text=recipe.text,
        date=str(recipe.date),
        json=recipe.json_data,
        code=recipe.template,
        image=encode_image(recipe.image) if show_image else '',
    )
    return store


def json_list(qs=None, show_image=False):
    output = {}
    projects = qs or Project.objects.all()
    for project in projects:
        proj_dict = encode_project(project, show_image=show_image)
        for recipe in project.analysis_set.all():
            proj_dict['recipes'][recipe.uid] = encode_recipe(recipe, show_image=show_image)
        output[project.uid] = proj_dict

    text = json.dumps(output, indent=4)
    return text


@api_error_wrapper(['GET'])
@ratelimit(key=RATELIMIT_KEY, rate='20/m')
def api_list(request):
    # Get the token and user
    token = auth.get_token(request=request)

    user = User.objects.filter(profile__token=token).first()

    # Get the project list corresponding to this user returns public projects if user is None.
    projects = auth.get_project_list(user=user)

    # Format the payload.
    payload = json_list(qs=projects)

    return HttpResponse(content=payload, content_type="text/json")


@api_error_wrapper(['GET', 'POST'])
#@token_access(klass=Project)
@csrf_exempt
@ratelimit(key='ip', rate='30/m')
def project_api(request, uid):
    """
    GET request : return project name, text, and image as a JSON file.
    POST request : change project name, text, and image given a JSON file.
    """

    qs = Project.objects.filter(uid=uid)

    # Get the json data with project info
    payload = json_list(qs, show_image=True)

    return HttpResponse(content=payload, content_type="text/json")


class Bunch(object):
    uid = ""
    name = ""
    text = ""
    date = ""
    privacy = ""
    image = ""
    recipes = []
    json = ""
    code = ""
    is_project = False

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def parse_json(json_dict):

    uid = json_dict.get('uid')
    name = json_dict.get('name')
    text = json_dict.get('text')
    date = json_dict.get('date')
    privacy = json_dict.get('privacy')
    image = json_dict.get('image')
    recipes = json_dict.get('recipes')
    json_text = json_dict.get('json')
    code = json_dict.get('code')
    is_project = recipes is not None
    data = Bunch(uid=uid, name=name, text=text, date=date, privacy=privacy,
                 image=image, recipes=recipes, is_project=is_project,
                 code=code, json=json_text)
    return data


def upload_recipe(obj, project, user=None, create=False):

    recipe = Analysis.objects.filter(uid=obj.uid).first()

    if not recipe and create:
        recipe = auth.create_analysis(project=project, user=user, uid=obj.uid,
                                      json_text=obj.json, template=obj.code)

    elif not recipe:
        return

    recipe.uid = obj.uid
    recipe.name = obj.name
    recipe.text = obj.text
    recipe.date = obj.date
    recipe.image = obj.image
    recipe.json_text = toml.dumps(obj.json)
    recipe.template = obj.code

    recipe.save()

    return


def upload_project(obj, user=None, create=False):

    # Check if this is a recipe or project

    project = Project.objects.filter(uid=obj.uid).first()

    if not project and create:
        project = auth.create_project(user=user, uid=obj.uid,
                                      name=obj.name,
                                      text=obj.text)
    elif not project:
        return

    project.uid = obj.uid
    project.name = obj.name
    project.text = obj.text
    project.date = obj.date
    project.privacy = obj.privacy
    project.image = obj.image

    project.save()
    for recipe, vals in obj.recipes.items():
        data = parse_json(vals)
        upload_recipe(data, project=project, user=user)

    return



@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Project)
@csrf_exempt
@ratelimit(key='ip', rate='30/m')
def api_upload(request):

    fstream = request.FILES.get("data", "")

    # Explicitly pass a 'create' flag
    create = request.POST.get('create', False)

    # Get a list of projects or a single one to update.
    json_obj = json.load(fstream)

    for key, item in json_obj.items():
        # Filter for uid if exists.
        data = parse_json(item)

        if data.is_project:
            upload_project(data, create=create)
        else:
            upload_recipe(data, create=create)

    return


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Data)
@csrf_exempt
@ratelimit(key='ip', rate='20/m')
def data_api(request):
    """
    GET request: Returns data
    PUT request: Updates file in data with given file.
    """

    uid = request.GET.get('uid', request.POST.get('uid'))
    data = Data.objects.filter(uid=uid).first()

    # Get the source that will replace target
    source = request.FILES.get("data", "")

    # Target first file in data directory.
    target = data.get_files()[0]

    if not target:
        msg = f"File does not exist."
        return HttpResponse(content=msg, content_type="text/plain")

    # Write source into target
    if request.method == "POST":
        # Validate source and target files before upload.
        valid, msg = auth.validate_file(source=source)
        if not valid:
            return HttpResponse(content=msg, content_type="text/plain")

        # Write source to target file.
        target = util.write_stream(stream=source, dest=target)

    # Return file contents in payload
    payload = open(target, 'r').read()

    return HttpResponse(content=payload, content_type="text/plain")


def job_api(request, uid):
    return
