
import hjson
import logging
import os

from django.conf import settings
from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from biostar.engine.models import Analysis, Project, image_path
from biostar.engine.decorators import require_api_key


logger = logging.getLogger("engine")


def get_thumbnail():
    return os.path.join(settings.STATIC_ROOT, "images", "placeholder.png")


def change_image(obj, file_object=None):

    if not obj:
        return get_thumbnail()

    obj.image.save(name=get_thumbnail(), content=file_object)

    return obj.image.path


@api_view(['GET'])
def project_api_list(request):

    projects = Project.objects.get_all()
    api_key = request.GET.get("k", "")

    # Only show public projects when api key is not correct or provided.
    if settings.API_KEY != api_key:
        projects = projects.filter(privacy=Project.PUBLIC)

    payload = []
    for project in projects:
        info = f"{project.uid}\t{project.name}\t{project.get_privacy_display()}\n"
        payload.append(info)

    payload = "".join(payload)
    return HttpResponse(content=payload, content_type="text/plain")


@api_view(['GET', 'PUT'])
@require_api_key(type=Project)
def project_info(request, uid):
    """
    GET request : return project info as json data
    PUT request : change project info using json data
    """

    project = Project.objects.get_all(uid=uid).first()

    if request.method == "PUT":
        file_object = request.data.get("file", "")
        conf = hjson.load(file_object)
        if file_object:
            project.name = conf.get("settings", {}).get("name") or project.name
            project.text = conf.get("settings", {}).get("help") or project.text
            project.save()

    payload = hjson.dumps(project.json_data, indent=4)

    return HttpResponse(content=payload, content_type="text/plain")


@api_view(['GET', 'PUT'])
@require_api_key(type=Project)
def project_image(request, uid):
    """
    GET request : return project image
    PUT request : change project image
    """
    project = Project.objects.filter(uid=uid).first()
    imgpath = project.image.path if project.image else get_thumbnail()

    if request.method == "PUT":
        file_object = request.data.get("file")
        imgpath = change_image(obj=project, file_object=file_object)

    data = open(imgpath, "rb") .read()

    return HttpResponse(content=data, content_type="image/jpeg")


@api_view(['GET', 'PUT'])
@require_api_key(type=Analysis)
def recipe_image(request, uid):
    """
    GET request: Return recipe image.
    PUT request: Updates recipe image with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()
    imgpath = recipe.image.path if recipe.image else get_thumbnail()

    if request.method == "PUT":
        file_object = request.data.get("file")
        imgpath = change_image(obj=recipe, file_object=file_object)

    data = open(imgpath, "rb").read()

    return HttpResponse(content=data, content_type="image/jpeg")


@api_view(['GET'])
def recipe_api_list(request, uid):

    api_key = request.GET.get("k", "")

    recipes = Analysis.objects.filter(project__uid=uid)
    # Only show public recipes when api key is not correct or provided.
    if settings.API_KEY != api_key:
        recipes = recipes.filter(project__privacy=Project.PUBLIC)

    payload = []
    for recipe in recipes:
        info = f"{recipe.uid}\t{recipe.name}\n"
        payload.append(info)

    payload = "".join(payload) if payload else "No recipes found."

    return HttpResponse(content=payload, content_type="text/plain")


@api_view(['GET', 'PUT'])
@require_api_key(type=Analysis)
def recipe_json(request, uid):
    """
    GET request: Returns recipe json
    PUT request: Updates recipe json with given file.
    """
    recipe = Analysis.objects.filter(uid=uid).first()

    if request.method == "PUT":
        # Get the new json that will replace the current one
        file_object = request.data.get("file", "")
        updated_json = hjson.load(file_object)
        recipe.json_text = hjson.dumps(updated_json) if file_object else recipe.json_text

        # Update help and name in recipe from json.
        if updated_json.get("settings"):
            recipe.name = updated_json["settings"].get("name", recipe.name)
            recipe.text = updated_json["settings"].get("help", recipe.text)

        recipe.save()

    payload = hjson.dumps(recipe.json_data, indent=4)

    return HttpResponse(content=payload, content_type="text/plain")


@api_view(['GET', 'PUT'])
@require_api_key(type=Analysis)
def recipe_template(request, uid):
    """
    GET request: Returns recipe template
    PUT request: Updates recipe template with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    # API key is always checked by @require_api_key decorator.
    if request.method == "PUT":
        # Get the new template that will replace the current one
        file_object = request.data.get("file", "")
        stream = file_object.read().decode("utf-8")
        recipe.template = stream if file_object else recipe.template
        recipe.save()
    payload = recipe.template

    return HttpResponse(content=payload, content_type="text/plain")

