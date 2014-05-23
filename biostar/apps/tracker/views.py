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

    seed_change = leech_change = completed_change = 0

    if event == 'started':
        if peer.left == 0:
            seed_change = 1
        else:
            leech_change = 1

    elif event == 'stopped':
        if peer.left == 0:
            seed_change = - 1
        else:
            leech_change = - 1
        peer.delete()

    elif event == 'completed':
        completed_change = 1

    try:
        logger.info('fetching %s' % info_hash)
        torrent = Torrent.objects.filter(info_hash=info_hash).update(
            uploaded=F('uploaded') + uploaded,
            downloaded=F('downloaded') + downloaded,
            seeds=F('seeds') + seed_change,
            leeches=F('leeches') + leech_change,
            completed=F('completed') + completed_change,
        )
    except Torrent.DoesNotExist:
        logger.inf('torrent does not exist for %s' % info_hash)

    # This below should should be cached

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


class TorrentUploadForm(forms.ModelForm):
    content = forms.FileField(label='File')

    class Meta:
        model = Torrent
        fields = ['name']


class TorrentUploadView(CreateView):
    model = Torrent
    template_name = 'tracker/torrent_upload.html'
    form_class = TorrentUploadForm
    success_url = reverse_lazy('tracker:index')

    def form_valid(self, form):
        form.instance.content = b''.join(form.cleaned_data['content'].chunks())
        return super(TorrentUploadView, self).form_valid(form)


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