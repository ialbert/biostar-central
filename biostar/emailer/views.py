from django.shortcuts import render
from .models import *
from django.db.models import Count

def index(request):

    groups = EmailGroup.objects.all().annotate(count=Count('subscription')).order_by("-count")
    emails = EmailAddress.objects.all().annotate(count=Count('subscription')).order_by("-count")[:20]

    subs = Subscription.objects.all()[:20]

    context = dict(groups=groups, emails=emails, subs=subs)
    return render(request, "emailer/index.html", context=context)
