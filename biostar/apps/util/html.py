import bleach

from django.conf import settings

from django.template import loader, Context, Template, RequestContext
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User
from re import match

# The pattern that matches the user link.
USER_PATTERN = "http://%s/u/\d+" % settings.SITE_DOMAIN

def userlinks(attrs, new=False):
    href = attrs['href']
    if match(USER_PATTERN, href):
        attrs['_text'] = "Test user"
    return attrs

# These callback will be applied on html parsing.
CALLBACKS = bleach.DEFAULT_CALLBACKS + [userlinks]

def render(name, **kwds):
    tmpl = loader.get_template(name)
    cont = Context(kwds)
    page = tmpl.render(cont)
    return page

