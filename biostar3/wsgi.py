"""
WSGI config for biostar3 project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.7/howto/deployment/wsgi/
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "biostar3.settings.base")

from django.core.wsgi import get_wsgi_application

def default():
    """
    Default WSGI application.
    """
    app = get_wsgi_application()
    return app

app = default()

def whitenoise():
    """
    Static content served via DjangoWhiteNoise.
    """
    from whitenoise.django import DjangoWhiteNoise
    app = get_wsgi_application()
    app = DjangoWhiteNoise(app)
    app.add_files('./export/media/', prefix='media/')
    return app


