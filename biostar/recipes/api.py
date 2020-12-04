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

    # Skip the placeholder images.
    #if text.startswith("iVBORw0KGgoAAAANSUhEUgAAAc8AAAHKCAIAAADq11fPAAAAAXNSR0IAr"):
    #    return ""

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
        image=encode_image(recipe.image) if show_image else ''
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
#@token_access(klass=Project, allow_create=True)
@csrf_exempt
@ratelimit(key='ip', rate='30/m')
def project_api(request, uid):
    """
    GET request : return project name, text, and image as a TOML file.
    POST request : change project name, text, and image given a TOML file.
    """

    qs = Project.objects.filter(uid=uid)

    token = auth.get_token(request=request)

    # Find the target user.
    user = User.objects.filter(profile__token=token).first()

    # Get the json data with project info
    payload = json_list(qs, show_image=True)

    #if request.method == "POST":
    #    # Fetch data from the
    #    stream = request.FILES.get("data")

    #    if stream:
    #        # Update or create a project using data.
    #        target = auth.update_project(obj=project, stream=stream,
    #                                     user=user, uid=uid,
    #                                     create=True, save=True)

    return HttpResponse(content=payload, content_type="text/json")


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Analysis, allow_create=True)
@csrf_exempt
@ratelimit(key='ip', rate='20/m')
def recipe_api(request):
    """
    GET request : return recipe json, template and image as a TOML string.
    POST request : change recipe json, template, and image given a TOML string.
    """

    # Get the object uid
    uid = request.GET.get('uid', request.POST.get('uid', ''))
    # Get the project uid in case of creation.
    pid = request.GET.get('pid', request.POST.get('pid', ''))

    recipe = Analysis.objects.filter(uid=uid).first()

    # Resolve the project from recipe or 'pid'
    project = recipe.project if recipe else None
    project = project or Project.objects.filter(uid=pid).first()

    target = recipe.api_data if recipe else {}

    if not project:
        return HttpResponse(content="Project does not exist.",
                            content_type="text/plain",
                            status=404)

    token = auth.get_token(request=request)
    # Find the target user.
    user = User.objects.filter(profile__token=token).first()

    # Replace source with target with valid POST request.
    if request.method == "POST":
        # Fetch data
        stream = request.FILES.get("data")
        if stream:
            # Update or create a recipe using data.
            target = auth.update_recipe(obj=recipe, stream=stream,
                                        save=True, create=True,
                                        user=user, uid=uid,
                                        project=project)
    # Get the payload as a toml file.
    payload = json.dumps(target)

    return HttpResponse(content=payload, content_type="text/plain")


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
