# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)

from django.db import models
from django.conf import settings


class Torrent(models.Model):
    name = models.CharField(max_length=200, default="Data", db_index=True)

    info_hash = models.CharField(max_length=40, db_index=True)

    # Some torrents may end up as disabled
    disabled = models.BooleanField(default=False)

    completed = models.PositiveIntegerField(default=0)
    downloaded = models.PositiveIntegerField(default=0)
    uploaded = models.PositiveIntegerField(default=0)

    # computed as len(Peer.objects.filter(torrent=self).filter(left=0)
    seeds = models.PositiveIntegerField(default=0)

    # computed as len(Peer.objects.filter(torrent=self).exclude(left=0)
    leeches = models.PositiveIntegerField(default=0)

    # Store the torrentin the database
    content = models.BinaryField()

    # Date related functionality
    lastupdate_date = models.DateTimeField(auto_now=True)
    creation_date = models.DateTimeField(auto_now_add=True)

    user = models.ForeignKey(settings.AUTH_USER_MODEL, null=True)

    class Meta:
        ordering = ('name', )

    def __unicode__(self):
        return self.name


class Peer(models.Model):
    """Contains peer fields"""

    torrent = models.ForeignKey(Torrent)

    peer_id = models.CharField(max_length=40, db_index=True)
    user_agent = models.CharField(max_length=120, blank=True)
    ip_address = models.IPAddressField()
    key = models.CharField(max_length=40, blank=True)
    port = models.IntegerField()
    uploaded = models.PositiveIntegerField(default=0)
    downloaded = models.PositiveIntegerField(default=0)
    left = models.PositiveIntegerField('remaining', default=0)

    lastupdate_date = models.DateTimeField(auto_now=True)
    creation_date = models.DateTimeField(auto_now_add=True)

    user = models.ForeignKey(settings.AUTH_USER_MODEL, null=True)

    def __unicode__(self):
        return "%s,%s" % (self.peerinfo, self.torrentinfo)


