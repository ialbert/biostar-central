# Create your views here.

from django import forms
from django.views.generic import  UpdateView
from django.contrib.flatpages.models import FlatPage
from django.http import HttpResponseRedirect
from django.contrib import messages

import logging

logger = logging.getLogger(__name__)

class FlatPageForm(forms.ModelForm):
    model = FlatPage
    content = forms.CharField(widget=forms.Textarea,
                              min_length=80, max_length=15000,
                              label="Post content")


class FlatPageUpdate(UpdateView):
    model = FlatPage
    fields = ['content']
    #form_class = FlatPageForm
    template_name = "flatpages/flatpage_edit.html"

    def post(self, *args, **kwargs):
        req = self.request
        user = req.user

        logger.info("user %s edited %s" % (user, kwargs))
        if not self.request.user.is_admin:
            logger.error("user %s access denied on %s" % (user, kwargs))
            messages.error(req, "Only administrators may edit that page")
            return HttpResponseRedirect("/")
        return super(FlatPageUpdate, self).post(*args, **kwargs)
