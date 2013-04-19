from django.views.generic import View, TemplateView, ListView, CreateView, UpdateView
from django.views.generic.base import TemplateResponseMixin
from django.core.urlresolvers import reverse
from django.shortcuts import render
from django.contrib import messages
from main.server import models
from django.db.models import Q
from main.server import html
from main.server.const import *
from django.conf import settings

class PageBase(TemplateView):
    url = "default"

    def get_context_data(self, **kwargs):
        context = super(PageBase, self).get_context_data(**kwargs)
        context['url'] = self.url
        return context

class MessageView(TemplateView):
    url = "/show/messages/"
    template_name = "refactored/message.page.html"

    def get_context_data(self, **kwargs):
        context = super(MessageView, self).get_context_data(**kwargs)

        user = self.request.user

        note_types = [ NOTE_USER, NOTE_PRIVATE ]

        q = Q(target=user, type__in=note_types)
        e = Q(sender=user, type=NOTE_USER)

        notes = models.Note.objects.filter(q).exclude(e)
        notes = notes.select_related('author', 'author__profile').order_by('-date')
        page = html.get_page(request=self.request, obj_list=notes, per_page=25)

        # reset the counts
        models.Note.objects.filter(target=user, unread=True).update(unread=False)
        models.UserProfile.objects.filter(user=user).update(new_messages=0)

        # the params object will carry
        layout  = settings.USER_PILL_BAR
        params  = html.Params(tab="", pill="messages", sort='', since='', layout=layout, title="Your Messages")

        context['page'] = page
        context['params'] = params

        return context

