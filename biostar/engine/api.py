
import hjson
import logging

from django.conf import settings
from django.http import HttpResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from biostar.engine.models import Analysis, Project
from biostar.utils.shortcuts import reverse
from biostar.engine.decorators import require_api_key


logger = logging.getLogger("engine")


@api_view(['GET'])
def project_api_list(request):

    projects = Project.objects.get_all()
    api_key = request.GET.get("k", "")

    # Only show public projects when api key is not correct or provided.
    if settings.API_KEY != api_key:
        projects = projects.filter(privacy=Project.PUBLIC)

    payload = dict()
    for project in projects:
        payload.setdefault(project.uid, dict()).update(
                            name=project.name,
                            recipes={recipe.uid:
                                     dict(name=recipe.name,
                                          json=reverse("recipe_api_json", kwargs=dict(uid=recipe.uid)),
                                          template=reverse("recipe_api_template", kwargs=dict(uid=recipe.uid)))
                                     for recipe in project.analysis_set.all()
                                     },
                            privacy=dict(Project.PRIVACY_CHOICES)[project.privacy],
                            )

    return Response(data=payload, status=status.HTTP_200_OK)


@api_view(['GET'])
def recipe_api_list(request):

    recipes = Analysis.objects.get_all()
    api_key = request.GET.get("k", "")

    # Only show public recipes when api key is not correct or provided.
    if settings.API_KEY != api_key:
        recipes = recipes.filter(project__privacy=Project.PUBLIC)

    payload = dict()
    for recipe in recipes:
        payload.setdefault(recipe.uid, dict()).update(
                            name=recipe.name,
                            json=reverse("recipe_api_json", kwargs=dict(uid=recipe.uid)),
                            template=reverse("recipe_api_template", kwargs=dict(uid=recipe.uid)),
                            privacy=dict(Project.PRIVACY_CHOICES)[recipe.project.privacy],
                            project_uid=recipe.project.uid,
                            project_name=recipe.project.name)

    return Response(data=payload, status=status.HTTP_200_OK)


@api_view(['GET', 'PUT'])
@require_api_key
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

    payload = recipe.json_data

    return Response(data=payload, status=status.HTTP_200_OK)


@api_view(['GET', 'PUT'])
@require_api_key
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


@api_view(['GET', 'PUT'])
@require_api_key
def recipe_image(request, uid):
    """
    GET request: Return recipe image.
    PUT request: Updates recipe image with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()
    img = recipe.image.path
    if request.method == "PUT":
        # Get the new image that will replace current one
        file_object = request.data.get("file", "")
        #TODO: check the file size?
        if file_object:
            stream = file_object.read()
            open(img, "wb").write(stream)

    img_stream = open(img, "rb").read()
    return HttpResponse(content=img_stream, content_type="image/jpeg")

