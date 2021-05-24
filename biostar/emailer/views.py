from django.shortcuts import render
from .models import *
from django.db.models import Count

def index(request):

    groups = EmailGroup.objects.all().annotate(count=Count('subscription')).order_by("-count")

    subs = EmailSubscription.objects.all()[:20]

    context = dict(groups=groups,subs=subs)
    return render(request, "emailer/index.html", context=context)
