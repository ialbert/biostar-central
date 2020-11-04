import toml
import json
import logging
import os
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


def tabular_list(qs=None):
    output = []
    projects = qs or Project.objects.all()
    for project in projects:
        for recipe in project.analysis_set.all():
            line = f"{project.uid}\t{project.name}\t{recipe.uid}\t{recipe.name}\t{project.get_privacy_display()}"
            output.append(line)

    return "\n".join(output)


@api_error_wrapper(['GET'])
@ratelimit(key=RATELIMIT_KEY, rate='20/m')
def api_list(request):

    # Get the token and user
    token = auth.get_token(request=request)
    user = User.objects.filter(profile__token=token).first()

    # Get the project list corresponding to this user
    # returns public projects if user is None.
    projects = auth.get_project_list(user=user)

    # Format the payload.
    payload = tabular_list(qs=projects)

    return HttpResponse(content=payload, content_type="text/plain")


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Project, allow_create=True)
@csrf_exempt
@ratelimit(key='ip', rate='20/m')
def project_api(request):
    """
    GET request : return project name, text, and image as a TOML file.
    POST request : change project name, text, and image given a TOML file.
    """
    # Get the object uid
    uid = request.GET.get('uid', request.POST.get('uid', ''))

    project = Project.objects.filter(uid=uid).first()

    token = auth.get_token(request=request)
    # Find the target user.
    user = User.objects.filter(profile__token=token).first()

    # Get the json data with project info
    target = project.api_data if project else {}

    if request.method == "POST":
        # Fetch data from the
        stream = request.FILES.get("data")

        if stream:
            # Update or create a project using data.
            target = auth.update_project(obj=project, stream=stream,
                                         user=user, uid=uid,
                                         create=True, save=True)

    payload = json.dumps(target)

    return HttpResponse(content=payload, content_type="text/plain")


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
        return HttpResponse(content="Project does not exist.", content_type="text/plain")

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

