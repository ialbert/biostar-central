
from .models import Analysis, Project
from django.conf import settings
from django.http import HttpResponse
from django.shortcuts import render


def recipe_list(request):

    recipes = Analysis.objects.get_all()
    api_key = request.GET.get("k", "")

    if settings.API_KEY == api_key:
        recipes = recipes.filter(project__privacy=Project.PUBLIC)

    context = dict(recipes=recipes)
    return render(request, "api_list.html", context=context)


def recipe_json(request, uid):
    """Returns json"""

    recipe = Analysis.objects.filter(uid=uid).first()
    api_key = request.GET.get("k", "")

    if settings.API_KEY == api_key:
        recipe = recipe.filter(project__privacy=Project.PUBLIC)

    return HttpResponse(recipe.json_text, content_type="application/json")


def recipe_template(request, uid):
    """Returns template"""

    recipe = Analysis.objects.filter(uid=uid).first()
    api_key = request.GET.get("k", "")

    if settings.API_KEY == api_key:
        recipe = recipe.filter(project__privacy=Project.PUBLIC)

    return HttpResponse(recipe.template, content_type="text/plain")


