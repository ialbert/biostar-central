# -*- encoding: utf-8 -*-
from __future__ import (unicode_literals, division, absolute_import)

from django.contrib import admin
from .models import *

admin.site.register(Torrent)
admin.site.register(Peer)

