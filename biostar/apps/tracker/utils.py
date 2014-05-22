# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)
from django.http import HttpResponse
from .bencode import bencode, bdecode
import hashlib

def torrent_get_length(content):
    decoded = bdecode(str(content))
    if 'length' in decoded['info']:
        return int(decoded['info']['length'])
    else:
        return sum(map(lambda x: x['length'], decoded['info']['files']))


def torrent_get_hash(content):
    decoded = bdecode(str(content))
    return hashlib.sha1(bencode(decoded['info'])).hexdigest()


def torrent_get_announce(content):
    decoded = bdecode(str(content))
    return str(decoded.get('announce'))


def torrent_set_announce(content, url):
    decoded = bdecode(str(content))
    decoded['announce'] = str(url)
    if 'announce-list' in decoded:
        del decoded['announce-list']
    return bencode(decoded)


def torrent_get_private(content):
    decoded = bdecode(str(content))
    return 'private' in decoded['info']


def torrent_set_private(content, private=True):
    decoded = bdecode(str(content))
    if private:
        decoded['info']['private'] = 1
    else:
        if 'private' in decoded['info']:
            del decoded['info']['private']
    return bencode(decoded)


def encode_peer_fast(peer):
    return 'd2:ip%d:%s4:porti%dee' % (len(peer.info.ip_address), peer.info.ip_address, peer.info.port)


def announce_response_fast(peers, complete, incomplete, interval=1800):
    return HttpResponse(
        'd8:intervali%de8:completei%de10:incompletei%de5:peersl%see' %
        (interval, complete, incomplete, "".join(map(encode_peer_fast, peers))),
        'text/plain')

