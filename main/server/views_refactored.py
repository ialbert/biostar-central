from django.views.generic import View, TemplateView, ListView, CreateView, UpdateView
from django.views.generic.base import TemplateResponseMixin
from django.core.urlresolvers import reverse
from django.shortcuts import render
from django.contrib import messages
from main.server import models
from django.db.models import Q, F
from main.server import html
from main.server.const import *
from django.conf import settings
from main import middleware
from datetime import datetime, timedelta
from django.shortcuts import render_to_response
from django.template import RequestContext
import random
from django.core.cache import cache

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

        note_types = [NOTE_USER, NOTE_PRIVATE]

        q = Q(target=user, type__in=note_types)
        e = Q(sender=user, type=NOTE_USER)

        notes = models.Note.objects.filter(q).exclude(e)
        notes = notes.select_related('author', 'author__profile').order_by('-date')
        page = html.get_page(request=self.request, obj_list=notes, per_page=25)

        # evaluate the query here so that the unread status is kept
        page.object_list = list(page.object_list)

        # reset the counts
        models.Note.objects.filter(target=user, unread=True).update(unread=False)
        models.UserProfile.objects.filter(user=user).update(new_messages=0)

        sess = middleware.Session(self.request)
        counts = sess.get_counts("message_count")
        sess.save()

        # the params object will carry
        layout = settings.USER_PILL_BAR

        params = html.Params(tab="", pill="messages",
                             sort='', since='', layout=layout, title="Your Messages")

        context['page'] = page
        context['params'] = params
        context['counts'] = counts

        return context

class AdHelp(TemplateView):
    url = "/help/ads/"
    template_name = "help/help.ads.html"

def get_ad_ids():
     ids = models.Ad.objects.filter(status=models.Ad.RUNNING).all().values_list("id")
     ids = [ x[0] for x in ids ]
     return ids

class NextAd(TemplateView):
    url = "/view/ad/"
    template_name = "refactored/ad.view.html"
    CACHE_KEY = "ads"

    def get_context_data(self, **kwargs):
        context = super(NextAd, self).get_context_data(**kwargs)

        ids = cache.get(self.CACHE_KEY)
        if not ids:
            ids = get_ad_ids()
            cache.set(self.CACHE_KEY, ids, 300)

        if ids:
            id = random.choice(ids)
            ad = models.Ad.objects.get(pk=id, status=models.Ad.RUNNING)
            ad.show_count += 1
            ad.save()
        else:
            ad = None

        context['ad'] = ad
        return context



def allow_start(user):

    if user.is_authenticated():

        ad_count = models.Ad.objects.filter(status=models.Ad.RUNNING, status_by=user).count()

        # minimize the database hits
        def func(ad, user):

            cond = user.profile.score > settings.AD_MIN_REP and ad_count < 1
            cond = cond or (user.profile.score > settings.AD_MOD_REP and ad_count < 2)
            return cond and (ad.status == models.Ad.STOPPED) and (ad.post.status == POST_OPEN)
    else:
        def func(ad, user):
            return False

    return func

def allow_stop(ad, user):
    if not user.is_authenticated():
        return False
    is_moderator = user.profile.can_moderate
    same_author = (user == ad.user)
    cond = (ad.status == models.Ad.RUNNING) and same_author

    return cond

class ToggleAd(ListView):
    url = "toggle-ad"

    def get(self, *args, **kwargs):

        url = reverse(AdView.url)

        m = models.Ad

        pk = kwargs['pk']
        user = self.request.user
        ad = models.Ad.objects.get(pk=pk)

        ad.start = allow_start(user)(ad=ad, user=user)
        ad.stop  = allow_stop(ad=ad, user=user)

        now = datetime.now()

        ad.expiration_date = now + timedelta(days=30)


        action = self.request.GET.get("action", "")
        if action == "start":
            if allow_start(user)(ad, user):
                ad.status_by = user
                ad.status = m.RUNNING
                ad.save()
                messages.info(self.request, 'Ad started.')
            else:
                messages.error(self.request, 'Ad starting was denied.')
            return html.redirect(url)

        if action == "stop":
            if allow_stop(ad, user):
                ad.status_by = user
                ad.status = m.STOPPED
                ad.save()
                messages.info(self.request, 'Ad stopped.')
            else:
                messages.error(self.request, 'Ad stopping was denied')
            return html.redirect(url)

        messages.error(self.request, 'Action not recognized')
        return html.redirect(url)


class AdView(ListView):
    url = "show-ads"
    template_name = "refactored/ad.list.html"
    paginate_by = 100
    context_object_name = 'ads'

    def get_search(self):
        term = self.request.GET.get("search", "")
        return term

    def get_queryset(self):

        user = self.request.user

        term = self.get_search()

        if term:
            cond = Q(user__profile__display_name__icontains=term) | Q(post__title__icontains=term) | \
                Q(status_by__profile__display_name__icontains=term)
            query = models.Ad.objects.filter(cond)
        else:
            query = models.Ad.objects.all()

        query = query.select_related("user__profile", "status_by__profile",
                                                                 "post").order_by('status', '-show_count')

        query = list(query[:30])
        for ad in query:
            ad.allow_start = allow_start(user)(ad=ad, user=user)
            ad.allow_stop = allow_stop(ad=ad, user=user)
            ad.ctr = float(ad.post.views)/(ad.show_count + 1) * 100

        return query

    def get_context_data(self, target="", **kwargs):
        context = super(AdView, self).get_context_data(**kwargs)
        user = self.request.user

        layout = settings.USER_PILL_BAR
        params = html.Params(tab="", pill="ads", sort='', since='', layout=layout, title="Ad List")

        context['params'] = params
        context['search'] = self.get_search()

        return context
