from django.shortcuts import render


def index(request):
    context = dict()
    return render(request, "emailer/index.html", context=context)
