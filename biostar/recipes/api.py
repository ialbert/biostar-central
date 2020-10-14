import toml as hjson
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


def tabular_list():
    output = []
    projects = Project.objects.all()
    for project in projects:
        for recipe in project.analysis_set.all():
            line = f"{project.uid}\t{project.name}\t{recipe.uid}\t{recipe.name}\t{project.get_privacy_display()}"
            output.append(line)

    return "\n".join(output)



@api_error_wrapper(['GET'])
@ratelimit(key=RATELIMIT_KEY, rate='20/m')
def api_list(request):
    payload = tabular_list()
    return HttpResponse(content=payload, content_type="text/plain")


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Project)
@ratelimit(key='ip', rate='20/m')
def project_api(request, uid):

    """
    GET request : return project name, text, and image as a TOML string.
    POST request : change project name, text, and image given a TOML string.
    """

    project = Project.objects.filter(uid=uid).first()

    # Get the json data with project info
    target = project.json_data

    # Replace source with target with valid POST request.
    if request.method == "POST":
        stream = request.FILES.get("file")
        source = hjson.load(stream.read())
        if stream:
            # Get new fields from the POST request and set them.
            target['text'] = project.text = source.get('text', project.text)
            target['name'] = project.name = source.get('name', project.name)
            target['help'] = project.help = source.get('help', project.help)

            #target['image'] = project.image = source.get('image', project.image)

            project.save()

    payload = hjson.dumps(target)

    return HttpResponse(content=payload, content_type="text/plain")


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Analysis)
@ratelimit(key='ip', rate='20/m')
def recipe_api(request, uid):
    """
    GET request : return recipe json, template and image as a TOML string.
    POST request : change recipe json, template, and image given a TOML string.
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    # Convert image to string
    img = auth.img_to_str(recipe.image) if recipe.image else auth.img_to_str(get_thumbnail())

    target = {"json": hjson.dumps(recipe.json_data),
              "template": recipe.template,
              "image": img}

    stream = request.FILES.get("file")
    source = hjson.load(stream.read())

    # Replace source with target with valid POST request.
    if request.method == "POST":
        # Get the toml object from the POST request
        interface = source.get('json', recipe.json_text)
        template = source.get('template', recipe.template)

        target['json'] = recipe.json_text = interface
        target['template'] = recipe.template = template


        # Swap the binary image
        #target['image'] = recipe.image = source.get('image', recipe.image)
        recipe.save()

    # Get the payload as a toml file.
    payload = hjson.dumps(target)

    return HttpResponse(content=payload, content_type="text/plain")


@api_error_wrapper(['GET', 'POST'])
@token_access(klass=Data)
@csrf_exempt
@ratelimit(key='ip', rate='20/m')
def data_api(request, uid):
    """
    GET request: Returns data
    PUT request: Updates file in data with given file.
    """

    data = Data.objects.filter(uid=uid).first()

    # Get the source that will replace target
    source = request.data.get("file", "")

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


#
# @api_error_wrapper(['GET', 'PUT'])
# @token_access(klass=Project)
# @ratelimit(key='ip', rate='20/m')
# def project_info(request, uid):
#     """
#     GET request : return project info as json data
#     PUT request : change project info using json data
#     """
#
#     project = Project.objects.filter(uid=uid).first()
#
#     if request.method == "PUT":
#         file_object = request.data.get("file", "")
#         conf = hjson.load(file_object)
#         if file_object:
#             project.name = conf.get("settings", {}).get("name") or project.name
#             project.text = conf.get("settings", {}).get("help") or project.text
#             project.save()
#
#     payload = json.dumps(project.json_data, indent=4)
#
#     return HttpResponse(content=payload, content_type="text/plain")
#

# @api_error_wrapper(['GET', 'PUT'])
# @token_access(klass=Project)
# @ratelimit(key='ip', rate='20/m')
# def project_image(request, uid):
#     """
#     GET request : return project image
#     PUT request : change project image
#     """
#     project = Project.objects.filter(uid=uid).first()
#     imgpath = project.image.path if project.image else get_thumbnail()
#
#     if request.method == "PUT":
#         file_object = request.data.get("file")
#         imgpath = change_image(obj=project, file_object=file_object)
#
#     data = open(imgpath, "rb").read()
#
#     return HttpResponse(content=data, content_type="image/jpeg")


# @api_error_wrapper(['GET', 'PUT'])
# @token_access(klass=Analysis)
# @ratelimit(key='ip', rate='20/m')
# def recipe_image(request, uid):
#     """
#     GET request: Return recipe image.
#     PUT request: Updates recipe image with given file.
#     """
#
#     recipe = Analysis.objects.filter(uid=uid).first()
#     imgpath = recipe.image.path if recipe.image else get_thumbnail()
#
#     if request.method == "PUT":
#         file_object = request.data.get("file")
#         imgpath = change_image(obj=recipe, file_object=file_object)
#
#     data = open(imgpath, "rb").read()
#
#     return HttpResponse(content=data, content_type="image/jpeg")

#
# @api_error_wrapper(['GET', 'PUT'])
# @token_access(klass=Analysis)
# @ratelimit(key='ip', rate='20/m')
# def recipe_json(request, uid):
#     """
#     GET request: Returns recipe json
#     PUT request: Updates recipe json with given file.
#     """
#     recipe = Analysis.objects.filter(uid=uid).first()
#
#     if request.method == "PUT":
#         # Get the new json that will replace the current one
#         file_object = request.data.get("file", "")
#         updated_json = hjson.load(file_object)
#         recipe.json_text = hjson.dumps(updated_json) if file_object else recipe.json_text
#
#         # Update help and name in recipe from json.
#         if updated_json.get("settings"):
#             recipe.name = updated_json["settings"].get("name", recipe.name)
#             recipe.text = updated_json["settings"].get("help", recipe.text)
#
#         recipe.save()
#
#     payload = json.dumps(recipe.json_data, indent=4)
#
#     return HttpResponse(content=payload, content_type="text/plain")


# @api_error_wrapper(['GET', 'PUT'])
# @token_access(klass=Analysis)
# @ratelimit(key='ip', rate='20/m')
# def recipe_template(request, uid):
#     """
#     GET request: Returns recipe template
#     PUT request: Updates recipe template with given file.
#     """
#
#     recipe = Analysis.objects.filter(uid=uid).first()
#
#     # API key is always checked by @require_api_key decorator.
#     if request.method == "PUT":
#         # Get the new template that will replace the current one
#         file_object = request.data.get("file", "")
#         stream = file_object.read().decode("utf-8")
#         recipe.template = stream if file_object else recipe.template
#         recipe.save()
#     payload = recipe.template
#
#     return HttpResponse(content=payload, content_type="text/plain")
#