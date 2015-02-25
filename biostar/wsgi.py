"""
WSGI config for biostar project.
"""

from django.core.wsgi import get_wsgi_application

# This is the default application
application = get_wsgi_application()

def white():
    # This is an alternative WSGI app that wraps static content
    from whitenoise.django import DjangoWhiteNoise
    white = get_wsgi_application()
    white = DjangoWhiteNoise(white)
    return white
