# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)

from django.db import models
from django.conf import settings
from django.contrib import admin

class Peer(models.Model):
    """Contains peer fields"""

    info_hash = models.CharField(max_length=40, db_index=True)
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


admin.site.register(Peer)
