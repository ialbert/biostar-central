# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)
from django.shortcuts import render, get_object_or_404
from django.db.models import F
from django.views.generic import ListView, DetailView, CreateView
from .models import Peer
from urlparse import parse_qs, urlsplit
from .bencode import bencode
from django.http import HttpResponse
import logging
from django.utils.text import slugify
from . import utils
from django.conf import settings
from django import forms
from django.core.urlresolvers import reverse_lazy
from biostar.apps.posts.models import Torrent
from django.views.generic.edit import DeleteView
from django.contrib import messages
from biostar.apps.tracker.utils import torrent_get_announce

logger = logging.getLogger(__name__)

# correct response bencode #
def bencoded_response(data):
    logger.info(data)
    return HttpResponse(bencode(data), content_type='text/plain')


def error_response(message):
    data = {b'failure reason': str(message)}
    return bencoded_response(data)


INTERVAL = settings.TRACKER_UPDATE_INTERVAL


def response(peers, complete, incomplete, interval=INTERVAL):
    encoded_peers = [{b'ip': str(peer.ip_address), b'port': int(peer.port)} for peer in peers]
    data = {
        b'interval': int(interval),
        b'complete': int(complete),
        b'incomplete': int(incomplete),
        b'peers': encoded_peers
    }
    return bencoded_response(data)


def scrape_response(torrentinfos):
    data = {b'files': torrentinfos}
    return bencoded_response(data)


def extract_params(request):
    """
    We can't use request.GET.get('info_hash') due django unicode craziness
    so we are parsing querystring by hand for 'info_hash' and 'peer_id'"""
    qs = urlsplit(request.get_full_path())[3]
    params = parse_qs(qs.encode('ascii'))
    info_hash = peer_id = ''
    if 'info_hash' in params and params['info_hash']:
        info_hash = params['info_hash'][0]
    if 'peer_id' in params and params['peer_id']:
        peer_id = params['peer_id'][0]
    return (info_hash, peer_id)


def create_torrent(request, info_hash):
    Torrent.objects.create(info_hash=info_hash.encode('hex'))

def create_magnet(request, pk):
    '''
    :returns: A magnet link conforming with magnet URI schema:
        magnet:?xt=urn:btih:'INFO_HASH'&dn='NAME'&tr='TRACKER1:port'&tr='TRACKER2:port'
    '''

    torrent = get_object_or_404(Torrent, pk=pk)

    magnet = "magnet:?xt=urn:btih:{info_hash}&dn={name}&tr={tracker}"
    magnet.format(info_hash=torrent.info_hash, name=torrent.name, tracker=torrent_get_announce(torrent.content))

    response = HttpResponse(magnet, content_type='text/plain')

    return response


def announce(request):
    """This function handles announce request send by torrent client.
    All input parametres are send as GET values.
    Function returns bencoded request with normal or error response in benc.
    """

    info_hash, peer_id = extract_params(request)

    # sanity_check info_hash
    if len(info_hash) != 20:
        return error_response('invalid info_hash')

    # sanity check peer_id
    if len(peer_id) != 20:
        return error_response('invalid peer_id')

    # sanity check port
    try:
        port = request.GET.get('port')
        port = int(port)
        if port < 1 or port > 65534:
            raise ValueError()
    except ValueError:
        return error_response('invalid port')

    # sanity check key
    key = request.GET.get('key', '')
    if len(key) > 100:
        return error_response('invalid key')


    # User throttling will be needed here
    info_hash = info_hash.encode('hex')


    # get peer info
    (peer, new_peer) = Peer.objects.get_or_create(
        info_hash=info_hash,
        ip_address=request.META.get('REMOTE_ADDR', '0.0.0.0'),
        port=port)

    # Update peer info if it was just created
    if new_peer:
        peer.peer_id = peer_id.encode('hex')
        peer.user_agent = request.META.get('HTTP_USER_AGENT', 'N/A'),
        peer.key = key,
        peer.save()

    # Get more parameters
    uploaded = int(request.GET.get('uploaded', 0))
    downloaded = int(request.GET.get('downloaded', 0))
    left = int(request.GET.get('left', 0))

    # Update stats for peers
    peer.uploaded = F('uploaded') + uploaded
    peer.downloaded = F('downloaded') + downloaded
    peer.left = left
    peer.save()

    # Handle the events
    event = request.GET.get('event', None)

    if event == 'stopped':
        peer.delete()

    try:
        logger.info('fetching %s' % info_hash)
        torrent = Torrent.objects.filter(info_hash=info_hash).update(
            uploaded=F('uploaded') + uploaded,
            seeds = Peer.objects.filter(info_hash=info_hash).count(),
        )

    except Torrent.DoesNotExist:
        logger.info('torrent does not exist for %s' % info_hash)

    # everything is ok, let's make response
    peers = Peer.objects.filter(info_hash=info_hash).order_by('?')[:20]

    # find the current stats on torrents
    query = Peer.objects.filter(info_hash=info_hash)
    complete = query.filter(left=0).count()
    incomplete = query.exclude(left=0).count()

    return response(peers, complete, incomplete)


def scrape(request):
    """This function handles scrape request of torrent which requests stats about torrent"""

    qs = urlsplit(request.get_full_path())[3]
    params = parse_qs(qs.encode('ascii'))
    info_hash_list = params.get('info_hash', [])

    if len(info_hash_list) == 0:
        return error_response('no info_hash specified, scraping all torrents is disabled')

    if len(info_hash_list) > 100:
        return error_response('too many info_hash specified')

    for info_hash in info_hash_list:
        if len(info_hash) != 20:
            return error_response('invalid info_hash')

    info_hash_list = map(lambda x: x.encode('hex'), info_hash_list)

    reply = {}
    for torrentinfo in Torrent.objects.filter(info_hash__in=info_hash_list):
        reply[torrentinfo.info_hash.decode('hex')] = {
            b'complete': torrentinfo.seeds,
            b'incomplete': torrentinfo.leeches,
            b'downloaded': torrentinfo.completed
        }

    return bencoded_response(reply)


def download_torrent(request, pk):
    torrent = get_object_or_404(Torrent, pk=pk)

    content = torrent.content

    #content = utils.torrent_set_announce(content, 'url')
    #content = utils.torrent_set_private(content, True)

    response = HttpResponse(content, content_type='application/octet-stream')
    response['Content-Disposition'] = 'attachment; filename="%s-%d.torrent"' % \
                                      (slugify(torrent.name), torrent.id)

    return response


class TorrentList(ListView):
    template_name = "tracker/torrent_list.html"
    queryset = Torrent.objects.select_related()
    paginate_by = 50
    context_object_name = "torrents"


def get_ip(request):
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or None
    return ip

class PeerList(ListView):
    template_name = "tracker/peer_list.html"
    queryset = Peer.objects.select_related()
    paginate_by = 50
    context_object_name = "peers"

    def get_context_data(self, **kwargs):
        context = super(PeerList, self).get_context_data(**kwargs)
        user = self.request.user
        ip = get_ip(self.request)
        if user.is_authenticated() and ip:
            peers = Peer.objects.filter(user=None, ip_address=ip).update(user=user)
            logger.info(peers)

        return context

class PeerDetail(DetailView):
    template_name = "tracker/peer_detail.html"
    model = Peer
    context_object_name = "peer"

class TorrentDetail(DetailView):
    template_name = "tracker/torrent_detail.html"
    model = Torrent
    context_object_name = "torrent"

    def get_context_data(self, **kwargs):

        from biostar.apps.posts.auth import post_permissions

        context = super(TorrentDetail, self).get_context_data(**kwargs)
        user = self.request.user
        ip = get_ip(self.request)
        if user.is_authenticated() and ip:
            peers = Peer.objects.filter(user=None, ip_address=ip).update(user=user)
            logger.info(peers)

        torrent = context[self.context_object_name]

        # find the post of the torrent
        post = torrent.post

        # adds the post permissions
        post = post_permissions(self.request, post)

        context['peers'] = Peer.objects.filter(info_hash=torrent.info_hash).all()[:100]
        context['post'] = post

        return context

from django.http import Http404

class TorrentDelete(DeleteView):
    model = Torrent

    def get_success_url(self):
        obj = super(TorrentDelete, self).get_object()
        return obj.post.get_absolute_url()

    def get_object(self, queryset=None):
        """ Hook to ensure object is owned by request.user. """
        from biostar.apps.posts.auth import post_permissions

        obj = super(TorrentDelete, self).get_object()
        post = post_permissions(self.request, obj.post)

        if not post.is_editable:
            messages.error(self.request, "Current user may not delete the tracker")
            raise Http404

        return obj


