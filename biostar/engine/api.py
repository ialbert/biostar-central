
import hjson

from django.conf import settings

from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response

from biostar.engine.models import Analysis, Project
from biostar.utils.shortcuts import reverse
from biostar.engine.decorators import require_api_key


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
                            json=reverse("api_json", kwargs=dict(uid=recipe.uid)),
                            template=reverse("api_template", kwargs=dict(uid=recipe.uid))
                            )

    return Response(data=payload, status=status.HTTP_200_OK)


@api_view(['GET', 'PUT'])
@require_api_key
def recipe_json(request, uid):
    """
    GET request: Returns recipe json
    PUT request: Updates recipe json with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()

    # API key is always checked by @require_api_key decorator.

    if request.method == "PUT":
        # Get the new json that will replace the current one
        file_object = request.data.get("file", "")
        recipe.json_text = hjson.dumps(hjson.load(file_object)) if file_object else recipe.json_text
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

    return Response(data=payload, status=status.HTTP_200_OK)


