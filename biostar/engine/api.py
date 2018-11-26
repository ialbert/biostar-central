
from .models import Analysis, Project
from django.conf import settings
from biostar.utils.shortcuts import reverse
from biostar.engine.decorators import require_api_key
import hjson
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response


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
    PUT request: Updates recipe template with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()
    if not recipe:
        return Response(data=dict(msg="Recipe does not exist."))

    # API is always required for put requests, checked by
    # @require_api_key decorator.
    if request.method == "PUT":
        # Get the new json that will replace the current one
        file_object = request.data.get("file", "")
        recipe.json_text = hjson.dumps(hjson.load(file_object)) if file_object else recipe.json_text
        recipe.save()
        payload = recipe.json_data

    else:
        # Only show public recipes when api key is not correct or provided.
        cond = recipe.project.is_private and settings.API_KEY == request.GET.get("k", "")
        msg = dict(msg="Private recipes can not be accessed without an API key param (?k=).")
        payload = recipe.json_data if (recipe.project.is_public or cond) else msg

    return Response(data=payload, status=status.HTTP_200_OK)


@api_view(['GET', 'PUT'])
@require_api_key
def recipe_template(request, uid):
    """
    GET request: Returns recipe template
    PUT request: Updates recipe template with given file.
    """

    recipe = Analysis.objects.filter(uid=uid).first()
    if not recipe:
        return Response(data=dict(msg="Recipe does not exist."))

    # API is always required for put requests, checked by
    # @require_api_key decorator.
    if request.method == "PUT":

        # Get the new template that will replace the current one
        file_object = request.data.get("file", "")
        stream = file_object.read().decode("utf-8")
        recipe.template = stream if file_object else recipe.template
        recipe.save()
        payload = dict(template=recipe.template)

    else:
        # Only show public recipes when api key is not correct or provided.
        cond = recipe.project.is_private and settings.API_KEY == request.GET.get("k", "")
        msg = dict(msg="Private recipes can not be accessed without an API key param (?k=).")
        payload = dict(template=recipe.template) if (recipe.project.is_public or cond) else msg

    return Response(data=payload, status=status.HTTP_200_OK)


