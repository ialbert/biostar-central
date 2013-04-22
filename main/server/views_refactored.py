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
from main import middleware
from datetime import datetime, timedelta


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


class ToggleAd(ListView):
    url = "/toggle/ad/"

    def get(self, *args, **kwargs):

        url = reverse(AdView.url, kwargs=dict(target="my"))

        m = models.Ad

        pk = kwargs['pk']
        user = self.request.user
        ad = models.Ad.objects.get(pk=pk)

        now = datetime.now()

        # these will only be saved on success
        ad.expiration_date = now + timedelta(days=30)
        ad.status_by = user

        # is this the author of the ad
        same_author = (user == ad.user)

        # does the author have sufficient reputation
        has_rep = (user.profile.score > settings.MINIMUM_AD_REP)

        # check the user moderation rights
        is_moderator = user.profile.can_moderate

        # how many active ads does the user have
        authored_ads = models.Ad.objects.filter(status=m.ACTIVE, user=user).count()
        approved_ads = models.Ad.objects.filter(status=m.ACTIVE, status_by=user).count()

        # how many ads does the user have
        total_ads = authored_ads + approved_ads

        # the condition for starting an ad
        user_start = has_rep or is_moderator

        # ad is active
        if ad.status == m.ACTIVE:

            if same_author:
                ad.status = m.STOPPED
                ad.save()
                messages.info(self.request, 'Ad stopped.')
                return html.redirect(url)

            messages.error(self.request, 'This active ad may not be changed by the current user')
            return html.redirect(url)

        # ad is pending
        if ad.status == m.PENDING:

            if total_ads:
                messages.error(self.request, 'Your aready have active ads. You must stop those before activating new ones')

            # a moderator may start at least one ad
            if is_moderator and not total_ads:
                ad.status = m.ACTIVE
                ad.status_by = user
                ad.save()
                messages.info(self.request, 'Ad activated.')
                return html.redirect(url)

            # the user has permission to start their ad
            if same_author and has_rep and not total_ads:
                ad.status = m.ACTIVE
                ad.status_by = user
                ad.save()
                messages.info(self.request, 'Ad activated.')
                return html.redirect(url)

            if same_author:
                ad.status = m.STOPPED
                ad.status_by = user
                ad.save()
                messages.info(self.request, 'Ad stopped. This ad will not be reviewed until placed on Pending status.')
                return html.redirect(url)

            messages.error(self.request, 'This pending ad may not be changed by the current user')
            return html.redirect(url)

        # the ad is stopped only the owner may put it into pending mode
        if ad.status == m.STOPPED:
            if same_author:
                ad.status = m.PENDING
                ad.status_by = user
                ad.save()
                messages.info(self.request, 'Ad placed in pending mode.')
                return html.redirect(url)

            messages.error(self.request, 'Only the author may put this ad in Pending mode.')
            return html.redirect(url)

        messages.error(self.request, 'You are not authorized to make changes for that ad. Please read the activation rules.')
        return html.redirect(url)


class AdView(ListView):
    url = "show/ads/"
    template_name = "refactored/show.ads.html"
    paginate_by = 25
    context_object_name = 'ads'

    def get_queryset(self):
        user = self.request.user
        target = self.kwargs['target']

        if target == "my":
            cond = Q(user=user) | Q(status_by=user)
            messages.info(self.request, 'Showing your ads. Switch to <a href="%s">all ads</a>' % reverse(self.url,
                                                                                                         kwargs=dict(
                                                                                                             target="all")))
        else:
            messages.info(self.request, 'Showing all ads. Switch to <a href="%s">your ads</a>' % reverse(self.url,
                                                                                                         kwargs=dict(
                                                                                                             target="my")))
            cond = Q()

        queryset = models.Ad.objects.filter(cond).select_related("user__profile", "status_by__profile",
                                                                 "post").order_by('-id')

        return queryset

    def get_context_data(self, target="", **kwargs):
        context = super(AdView, self).get_context_data(**kwargs)
        user = self.request.user

        layout = settings.USER_PILL_BAR
        params = html.Params(tab="", pill="ads", sort='', since='', layout=layout, title="Ad List")

        sess = middleware.Session(self.request)
        counts = sess.get_counts("ad_count")
        sess.save()

        context['params'] = params
        context['user'] = self.request.user
        context['counts'] = counts
        context['ads'] = context['ads']

        return context
