import bleach

from django.conf import settings
from django.template import loader, Context, Template, RequestContext
import re
import bleach
import logging

logger = logging.getLogger(__name__)

ALLOWED_TAGS = bleach.ALLOWED_TAGS + settings.ALLOWED_TAGS
ALLOWED_STYLES = bleach.ALLOWED_STYLES + settings.ALLOWED_STYLES
ALLOWED_ATTRIBUTES = dict(bleach.ALLOWED_ATTRIBUTES)
ALLOWED_ATTRIBUTES.update(settings.ALLOWED_ATTRIBUTES)

# The pattern that matches the user link.
USER_PATTERN = "http://%s/u/(?P<uid>(\d+))" % settings.SITE_DOMAIN
POST_PATTERN = "http://%s/p/(?P<uid>(\d+))" % settings.SITE_DOMAIN

# These are for debugging only.
#USER_PATTERN = "http://%s/u/(?P<uid>(\d+))" % "abc.com"
#POST_PATTERN = "http://%s/p/(?P<uid>(\d+))" % "abc.com"

USER_RE = re.compile(USER_PATTERN)
POST_RE = re.compile(POST_PATTERN)

def parse_html(text):
    from biostar.apps.users.models import User
    from biostar.apps.posts.models import Post

    def postlinks(attrs, new=False):
        "Matches a user"
        try:
            href = attrs['href']
            patt = POST_RE.search(href)
            if patt:
                uid = patt.group("uid")
                attrs['_text'] = Post.objects.get(id=uid).title
        except Exception, exc:
            logger.error(exc)
        return attrs

    def userlinks(attrs, new=False):
        "Matches a user"
        try:
            href = attrs['href']
            patt = USER_RE.search(href)
            if patt:
                uid = patt.group("uid")
                attrs['_text'] = User.objects.get(id=uid).name
        except Exception, exc:
            logger.error(exc)
        return attrs

    CALLBACKS = bleach.DEFAULT_CALLBACKS + [userlinks, postlinks]

    text = bleach.linkify(text, callbacks=CALLBACKS, skip_pre=True)

    html = bleach.clean(text, tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)

    return html

def strip_tags(text):
    "Strip html tags from an input"
    text = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
    return text

def render(name, **kwds):
    tmpl = loader.get_template(name)
    cont = Context(kwds)
    page = tmpl.render(cont)
    return page

def test():
    pass

if __name__ == '__main__':
    test()
