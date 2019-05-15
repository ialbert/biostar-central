from django.shortcuts import render
from .models import *

def index(request):
    groups = EmailGroup.objects.all()[:20]
    emails = EmailAddress.objects.all().select_related("group")[:20]
    context = dict(groups=groups, emails=emails)
    return render(request, "emailer/index.html", context=context)
