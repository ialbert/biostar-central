from __future__ import absolute_import

import os
from celery import Celery
from django.conf import settings

app = Celery('biostar3')

app.config_from_object('django.conf:settings')

app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)

