import hjson
from .models import Analysis, Project
from django.conf import settings
from django.http import HttpResponse
from biostar.utils.shortcuts import reverse


def recipe_list(request):

    recipes = Analysis.objects.get_all()
    api_key = request.GET.get("k", "")

    # Only show public recipes when api key is not correct or provided.
    if settings.API_KEY != api_key:
        recipes = recipes.filter(project__privacy=Project.PUBLIC)

    payload = dict()
    for recipe in recipes:

        data = dict(name=recipe.name,
                    json=reverse("api_json", kwargs=dict(uid=recipe.uid)),
                    template=reverse("api_template", kwargs=dict(uid=recipe.uid))
                    )
        payload.setdefault(recipe.uid, dict()).update(data)

    payload = hjson.dumps(payload)
    response = HttpResponse(payload, content_type="application/json")

    return response


def recipe_json(request, uid):
    """
    Returns json
    """

    recipe = Analysis.objects.filter(uid=uid)
    api_key = request.GET.get("k", "")

    # Only show public recipes when api key is not correct or provided.
    if settings.API_KEY != api_key:
        recipe = recipe.filter(project__privacy=Project.PUBLIC)

    data = recipe.first().json_text if recipe else "Recipe does not exist."

    return HttpResponse(data, content_type="application/json")


def recipe_template(request, uid):
    """
    Returns template
    """

    recipe = Analysis.objects.filter(uid=uid)
    api_key = request.GET.get("k", "")

    # Only show public recipes when api key is not correct or provided.
    if settings.API_KEY != api_key:
        recipe = recipe.filter(project__privacy=Project.PUBLIC)

    data = recipe.first().template if recipe else "Recipe does not exist."
    return HttpResponse(data, content_type="text/plain")


