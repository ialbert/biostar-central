__author__ = 'ialbert'

from django.views.generic import View, TemplateView, ListView, CreateView, UpdateView
from django.views.generic.base import TemplateResponseMixin
from django.core.urlresolvers import reverse
from django.shortcuts import render
from django.contrib import messages
from main.server import models
from django.db.models import Q, F
from django.http import HttpResponseRedirect

def replace_tags(request, pairs):
    """Replaces a tag with a new value"""

    for id, value in pairs:
        value = value.strip()

        if not value:
            continue

        tag1 = models.Tag.objects.get(id=id)
        tag2, flag = models.Tag.objects.get_or_create(name=value)

        if tag1.name == tag2.name:
            continue
        if tag1.count > 25:
            messages.info(request, "Cannot replace tag: %s It has more than 25 uses" % tag1.name)
            return

        msg = "Replaced tag %s (%s) with %s (%s)" % (tag1.name, tag1.count, tag2.name, tag2.count)
        messages.info(request, msg)
        posts = tag1.post_set.all()
        for post in posts:
            post.tag_val = post.tag_val.replace(tag1.name, tag2.name)
            post.set_tags()

        tag1.delete()

def delete_tags(request, ids):

    tags = models.Tag.objects.filter(id__in=ids)
    for tag in tags:
        if tag.count > 5:
            messages.warning(request, "Tag %s has been used %s times - it should probably be merged/renamed" % (tag.name, tag.count))
            continue
        posts = tag.post_set.all()
        for post in posts:
            post.tag_val = post.tag_val.replace(tag.name, "")
            post.set_tags()
        tag.delete()

    messages.info(request, "Deleted tags: %s" % ', '.join(t.name for t in tags))

class TagList(TemplateView):
    url = "list-tags"
    template_name = "refactored/full.tag.list.html"

    def post(self,  *args, **kwargs):
        user = self.request.user
        is_moderator = not user.is_anonymous() and user.profile.can_moderate

        delete_ids = self.request.POST.getlist("delete")
        replace_pairs = [ (key[4:], value) for (key, value) in self.request.POST.items() if key.startswith("rep-") ]

        if is_moderator and replace_pairs:
            replace_tags(self.request, replace_pairs)

        if is_moderator and delete_ids:
            delete_tags(self.request, delete_ids)

        return HttpResponseRedirect(reverse(self.url))

    def get(self,  *args, **kwargs):
        user = self.request.user
        tags = self.get_queryset()
        is_moderator = not user.is_anonymous() and user.profile.can_moderate
        if not is_moderator:
            messages.warning(self.request, "Sorry you must be a moderator to edit tags")
            return HttpResponseRedirect("/")
        return render(self.request, self.template_name, {'tags': tags})

    def get_queryset(self):
        return models.Tag.objects.filter().order_by("name")