
from django.http import JsonResponse
from .models import  Analysis



def recipe_details(request, uid):
    "Returns recipe details that can be used to build another one."

    recipe = Analysis.objects.filter(uid=uid).first()

    data = dict(uid=recipe.uid, deleted=recipe.deleted,
                sticky=recipe.sticky, name=recipe.name,
                summary=recipe.summary, text=recipe.text,
                owner=recipe.owner.email, diff_author=recipe.diff_author.email,
                diff_date=recipe.diff_date, json_text=recipe.json_text,
                template=recipe.template, last_valid=recipe.last_valid,
                date=recipe.date)

    return JsonResponse(data)