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

# Matching patterns will be filled in with post title or user name
USER_PATTERN = r"http(s)?://%s/u/(?P<uid>(\d+))" % settings.SITE_DOMAIN
POST_PATTERN1 = r"http(s)?://%s/p/(?P<uid>(\d+))" % settings.SITE_DOMAIN
POST_PATTERN2 = r"http(s)?://%s/p/\d+/\#(?P<uid>(\d+))" % settings.SITE_DOMAIN

# Matches gists that may be embeded
GIST_PATTERN = r"https://gist.github.com/(?P<uid>([\w/]+))"

# Matches Youtube video links.
YOUTUBE_PATTERN = r"http(s)?://www.youtube.com/watch\?v=(?P<uid>(\w+))"

USER_RE = re.compile(USER_PATTERN)
POST_RE1 = re.compile(POST_PATTERN1)
POST_RE2 = re.compile(POST_PATTERN2)
GIST_RE = re.compile(GIST_PATTERN)
YOUTUBE_RE = re.compile(YOUTUBE_PATTERN)

def clean(text):
    "Sanitize text with no other substitutions"
    html = bleach.clean(text, tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)
    return html

def parse_html(text):
    "Sanitize text and expand links to match content"
    from biostar.apps.users.models import User
    from biostar.apps.posts.models import Post

    # This will collect the objects that could be embedded
    embed = []

    def internal_links(attrs, new=False):
        "Matches a user"
        try:
            href = attrs['href']

            # Try the patterns
            patt1 = POST_RE1.search(href)
            patt2 = POST_RE2.search(href)
            patt = patt1 or patt2
            if patt:
                uid = patt.group("uid")
                attrs['_text'] = Post.objects.get(id=uid).title

            # Try the user patterns
            patt3 = USER_RE.search(href)
            if patt3:
                uid = patt3.group("uid")
                attrs['_text'] = User.objects.get(id=uid).name

        except Exception, exc:
            logger.error(exc)
        return attrs

    def embedder(attrs, new):
        # This is an existing <a> tag, leave it be.
        if not new:
            return attrs

        href = attrs['_text']

        # Don't linkify non http links
        if href[:4] not in ('http', 'ftp:'):
            return None

        # Try the gist embedding patterns
        targets = [
            (GIST_RE, '<script src="https://gist.github.com/%s.js"></script>'),
            (YOUTUBE_RE, '<iframe width="420" height="315" src="//www.youtube.com/embed/%s" frameborder="0" allowfullscreen></iframe>')
        ]

        for regex, text in targets:
            patt = regex.search(href)
            if patt:
                uid = patt.group("uid")
                obj = text % uid
                embed.append( (uid, obj) )
                attrs['_text'] = uid
                attrs['href'] = uid
                if 'rel' in attrs:
                    del attrs['rel']

        return attrs

    CALLBACKS = bleach.DEFAULT_CALLBACKS + [embedder, internal_links]

    html = bleach.clean(text, tags=ALLOWED_TAGS,
        attributes=ALLOWED_ATTRIBUTES, styles=ALLOWED_STYLES)

    try:
        html = bleach.linkify(html, callbacks=CALLBACKS, skip_pre=True)
        # embed the objects
        for uid, obj in embed:
            emb_patt = '<a href="%s">%s</a>' % (uid, uid)
            html = html.replace(emb_patt, obj)
    except Exception, exc:
        logger.error("*** %s" % exc)

    return html

def strip_tags(text):
    "Strip html tags from text"
    text = bleach.clean(text, tags=[], attributes=[], styles={}, strip=True)
    return text

def render(name, **kwds):
    "Helper function to render a template"
    tmpl = loader.get_template(name)
    cont = Context(kwds)
    page = tmpl.render(cont)
    return page

def test():
    pass

if __name__ == '__main__':
    test()
