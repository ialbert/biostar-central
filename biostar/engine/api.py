
from .models import Analysis
from django.http import HttpResponse
from .const import *


def api_text(recipe, text_type):

    action_mapper = {JSON: recipe.json_text, TEMPLATE: recipe.template}
    # Only show public recipes for now.
    text = action_mapper.get(text_type, "") if recipe and recipe.project.is_public else ""

    return text


def api_recipe_json(request, uid):
    """Returns plain text version of json"""

    recipe = Analysis.objects.filter(uid=uid).first()
    text = api_text(recipe=recipe, text_type=JSON)

    return HttpResponse(text, content_type="text/plain")


def api_recipe_template(request, uid):
    """Returns plain text version of template"""

    recipe = Analysis.objects.filter(uid=uid).first()
    text = api_text(recipe=recipe, text_type=TEMPLATE)

    return HttpResponse(text, content_type="text/plain")


