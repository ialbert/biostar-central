"""
Markdown parser to render the Biostar style markdown.
"""
import re
import inspect, logging
from functools import partial
import mistune
import requests
from xml.sax.saxutils import unescape
from django.shortcuts import reverse
from django.db.models import F
import bleach
from bleach.linkifier import Linker
from bleach.linkifier import LinkifyFilter
import html2text
from django.conf import settings
from mistune import Renderer, InlineLexer, InlineGrammar
from mistune import escape as escape_text
from bleach.sanitizer import Cleaner
from biostar.forum import auth
from biostar.forum.models import Post, Subscription
from biostar.accounts.models import Profile, User
from bleach.callbacks import nofollow

logger = logging.getLogger('engine')

# Test input.
TEST_INPUT = '''

https://www.biostars.org/p/1 https://www.biostars.org/p/2

http://localhost:8000/p/3/

http://localhost:8000/accounts/profile/5

@test

<b>BOLD</b>

https://www.youtube.com/watch?v=Hc8QdwfYFT8

```
print 123
http://www.psu.edu
```

https://twitter.com/Linux/status/2311234267

'''

TEST_INPUT2 = '''


http://test.biostars.org/accounts/profile/user-2/ 

http://test.biostars.org/p/p371285/

'''

# Shortcut to re.compile
rec = re.compile

# Biostar patterns
PORT = ':' + settings.HTTP_PORT if settings.HTTP_PORT else ''
SITE_URL = f"{settings.SITE_DOMAIN}{PORT}"
USER_PATTERN = rec(fr"http(s)?://{SITE_URL}/u/(?P<uid>[\w_.-]+)(/)?")
POST_TOPLEVEL = rec(fr"http(s)?://{SITE_URL}/p/(?P<uid>(\w+))(/)?")
POST_ANCHOR = rec(fr"http(s)?://{SITE_URL}/p/\w+/\#(?P<uid>(\w+))(/)?")

# Match any alphanumeric characters after the @.
# These characters are allowed in handles: _  .  -
MENTINONED_USERS = rec(r"(\@(?P<handle>[\w_.'-]+))")

ALLOWED_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],
}

ALLOWED_TAGS = ['p', 'div', 'br', 'code', 'pre', 'h1', 'h2', 'h3', 'h4', 'hr', 'span', 's',
                'sub', 'sup', 'b', 'i', 'img', 'strong', 'strike', 'em', 'underline',
                'super', 'table', 'thead', 'tr', 'th', 'td', 'tbody', 'del', 'details', 'summary']
ALLOWED_TAGS = bleach.ALLOWED_TAGS + ALLOWED_TAGS

ALLOWED_PROTOCOLS = ['ftp', 'http']
ALLOWED_PROTOCOLS = bleach.ALLOWED_PROTOCOLS + ALLOWED_PROTOCOLS

ALLOWED_STYLES = ['color', 'font-weight', 'background-color', 'width height']

# Youtube patterns
# https://www.youtube.com/watch?v=G7RDn8Xtf_Y
YOUTUBE_PATTERN1 = rec(r"^http(s)?://www.youtube.com/watch\?v=(?P<uid>([\w-]+))(/)?")
YOUTUBE_PATTERN2 = rec(r"^https://www.youtube.com/embed/(?P<uid>([\w-]+))(/)?")
YOUTUBE_PATTERN3 = rec(r"^https://youtu.be/(?P<uid>([\w-]+))(/)?")
YOUTUBE_HTML = '<iframe width="420" height="315" src="//www.youtube.com/embed/%s" frameborder="0" allowfullscreen></iframe>'

# Ftp link pattern.
FTP_PATTERN = rec(r"^ftp://[\w\.]+(/?)$")

# Gist pattern.
# https://gist.github.com/afrendeiro/6732a46b949e864d6803
GIST_PATTERN = rec(r"^https://gist.github.com/(?P<uid>([\w/]+))")
GIST_HTML = '<script src="https://gist.github.com/%s.js"></script>'

# Twitter pattern.
# https://twitter.com/Linux/status/2311234267
TWITTER_PATTERN = rec(r"http(s)?://(www)?.?twitter.com/\w+/status(es)?/(?P<uid>([\d]+))(/)?([^\s]+)?")


def get_tweet(tweet_id):
    """
    Get the HTML code with the embedded tweet.
    It requires an API call at https://api.twitter.com/1/statuses/oembed.json as documented here:
    https://dev.twitter.com/docs/embedded-tweets - section "Embedded Tweets for Developers"
    https://dev.twitter.com/docs/api/1/get/statuses/oembed
    Params:
    tweet_id -- a tweet's numeric id like 2311234267 for the tweet at
    https://twitter.com/Linux/status/2311234267
    """
    try:
        response = requests.get("https://api.twitter.com/1/statuses/oembed.json?id={}".format(
            tweet_id))
        return response.json()['html']
    except:
        return ''


class MonkeyPatch(InlineLexer):
    """
    Mistune uses class attributes for default_rules and those do provide isolation
    between different instances of the parser.
    This subclass moves the default_rules to instance attributes.
    """

    def __init__(self, *args, **kwds):
        super(MonkeyPatch, self).__init__(*args, **kwds)
        self.default_rules = list(InlineLexer.default_rules)


class BiostarInlineGrammer(InlineGrammar):
    """
    Mistune has a rule 'text' that sees '@' symbol as part of a regular string.
    This excludes characters around the '@' symbol from being seen by other rules.
    This subclass overrides the 'text' rule to contain `@` as a stop element in the lookahead assertion

    Also the '_' character not considered special use * instead
    """
    text = re.compile(r'^[\s\S]+?(?=[\\<!\[*`~@]|https?://| {2,}\n|$)')


class BiostarRenderer(Renderer):

    def codespan(self, text):
        """Rendering inline `code` text.

        :param text: text content for inline code.
        """
        text = escape_text(text.rstrip(), smart_amp=True)
        return '<code>%s</code>' % text

    def block_code(self, code, lang=None):
        """
        This is overrides to turn smart_amp=True
        Rendering block level code. ``pre > code``.

        :param code: text content of the code block.
        :param lang: language of the given code.

        Turn smart_amp=True here to prevent &gt; changing to &amp;gt; after bleach clean.

        """
        code = code.rstrip('\n')
        if not lang:
            code = escape_text(code, smart_amp=True)
            return '<pre><code>%s\n</code></pre>\n' % code
        code = escape_text(code, quote=True, smart_amp=True)
        return '<pre><code class="lang-%s">%s\n</code></pre>\n' % (lang, code)


def rewrite_static(link):
    # Link is already a full path or external
    if link.startswith("/") or link.startswith("http"):
        return link

    # Make the link absolute to the static url
    link = "/static/" + link

    return link


class BiostarInlineLexer(MonkeyPatch):
    grammar_class = BiostarInlineGrammer

    def __init__(self, root=None, allow_rewrite=False, *args, **kwargs):
        """
        :param root: Root post that is being pared
        :param static_imgs:
        """
        self.root = root
        self.allow_rewrite = allow_rewrite

        super(BiostarInlineLexer, self).__init__(*args, **kwargs)
        self.enable_all()

    def enable_all(self):
        self.enable_post_link()
        self.enable_mention_link()
        self.enable_anchor_link()
        self.enable_user_link()
        self.enable_gist_link()
        self.enable_youtube_link1()
        self.enable_youtube_link2()
        self.enable_youtube_link3()
        self.enable_ftp_link()
        self.enable_twitter_link()

    def enable_post_link(self):
        self.rules.post_link = POST_TOPLEVEL
        self.default_rules.insert(0, 'post_link')

    def enable_mention_link(self):
        self.rules.mention_link = MENTINONED_USERS
        self.default_rules.insert(0, 'mention_link')

    def _process_link(self, m, link, title=None):
        line = m.group(0)
        text = m.group(1)
        if line[0] == '!':
            if self.allow_rewrite:
                # Ensure the link is a full url path found in to static directory.
                link = rewrite_static(link)

            return self.renderer.image(link, title, text)

        self._in_link = True
        text = self.output(text)
        self._in_link = False
        return self.renderer.link(link, title, text)

    def output_mention_link(self, m):

        handle = m.group("handle")
        # Query user and get the link
        user = User.objects.filter(profile__handle=handle).first()
        if user:
            profile = reverse("user_profile", kwargs=dict(uid=user.profile.uid))
            link = f'<a href="{profile}">{user.profile.name}</a>'
            # Subscribe mentioned users to post.
            if self.root:
                # Create user subscription if it does not already exist.
                auth.create_subscription(post=self.root, user=user, update=True)
        else:
            link = m.group(0)

        return link

    def output_post_link(self, m):
        uid = m.group("uid")
        link = m.group(0)
        post = Post.objects.filter(uid=uid).first()
        title = post.root.title if post else "Post not found"
        return f'<a href="{link}">{title}</a>'

    def enable_anchor_link(self):
        self.rules.anchor_link = POST_ANCHOR
        self.default_rules.insert(0, 'anchor_link')

    def output_anchor_link(self, m):
        uid = m.group("uid")
        link = m.group(0)
        post = Post.objects.filter(uid=uid).first()
        title = post.root.title if post else "Post not found"
        return f'<a href="{link}">{title}</a>'

    def enable_user_link(self):
        self.rules.user_link = USER_PATTERN
        self.default_rules.insert(0, 'user_link')

    def output_user_link(self, m):
        uid = m.group("uid")
        link = m.group(0)
        profile = Profile.objects.filter(uid=uid).first()
        name = profile.name if profile else f"Invalid user uid: {uid}"
        return f'<a href="{link}">{name}</a>'

    def enable_youtube_link1(self):
        self.rules.youtube_link1 = YOUTUBE_PATTERN1
        self.default_rules.insert(1, 'youtube_link1')

    def output_youtube_link1(self, m):
        uid = m.group("uid")
        # Isolate links to be later dealt with
        return m.group()

    def enable_youtube_link2(self):
        self.rules.youtube_link2 = YOUTUBE_PATTERN2
        self.default_rules.insert(1, 'youtube_link2')

    def output_youtube_link2(self, m):
        uid = m.group("uid")
        # Isolate links to be later dealt with
        return m.group()

    def enable_youtube_link3(self):
        self.rules.youtube_link3 = YOUTUBE_PATTERN3
        self.default_rules.insert(1, 'youtube_link3')

    def output_youtube_link3(self, m):
        uid = m.group("uid")
        # Isolate links to be later dealt with
        return m.group()

    def enable_twitter_link(self):
        self.rules.twitter_link = TWITTER_PATTERN
        self.default_rules.insert(1, 'twitter_link')

    def output_twitter_link(self, m):
        uid = m.group("uid")
        # Isolate links to be later dealt with
        return m.group()

    def enable_gist_link(self):
        self.rules.gist_link = GIST_PATTERN
        self.default_rules.insert(3, 'gist_link')

    def output_gist_link(self, m):
        uid = m.group("uid")
        # Isolate links to be later dealt with
        return m.group()

    def enable_ftp_link(self):
        self.rules.ftp_link = FTP_PATTERN
        self.default_rules.insert(3, 'ftp_link')

    def output_ftp_link(self, m):
        link = m.group(0)
        return f'<a href="{link}">{link}</a>'


def embedder(attrs, new, embed=None):
    embed = [] if embed is None else embed

    # Existing <a> tag, leave as is.
    if not new:
        return attrs

    href = attrs['_text']
    linkable = href.startswith(('http', 'ftp:', 'https'))

    # Don't linkify non http links
    if not linkable:
        return None

    # Try embedding patterns
    targets = [
        (GIST_PATTERN, lambda x: GIST_HTML % x),
        (YOUTUBE_PATTERN1, lambda x: YOUTUBE_HTML % x),
        (YOUTUBE_PATTERN2, lambda x: YOUTUBE_HTML % x),
        (YOUTUBE_PATTERN3, lambda x: YOUTUBE_HTML % x),
        (TWITTER_PATTERN, get_tweet),
    ]

    for regex, get_text in targets:
        patt = regex.search(href)
        if patt:
            uid = patt.group("uid")
            obj = get_text(uid)
            embed.append((patt.group(), obj))
            attrs['_text'] = patt.group()
            if 'rel' in attrs:
                del attrs['rel']

    return attrs


def linkify(text):
    # List of links to embed
    embed = []
    html = bleach.linkify(text=text, callbacks=[partial(embedder, embed=embed), nofollow], skip_tags=['pre', 'code'])

    # Embed links into html.
    for em in embed:
        source, target = em
        emb = f'<a href="{source}" rel="nofollow">{source}</a>'
        html = html.replace(emb, target)

    return html


def safe(f):
    """
    Safely call an object without causing errors
    """
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as exc:
            logger.error(f"Error with {f.__name__}: {exc}")
            text = kwargs.get('text', args[0])
            return text

    return inner


@safe
def parse(text, post=None, clean=True, escape=True, allow_rewrite=False):
    """
    Parses markdown into html.
    Expands certain patterns into HTML.

    clean : Applies bleach clean BEFORE mistune escapes unsafe characters.
            Also removes unbalanced tags at this stage.
    escape  : Escape html originally found in the markdown text.
    allow_rewrite : Serve images with relative url paths from the static directory.
                  eg. images/foo.png -> /static/images/foo.png
    """

    # Resolve the root if exists.
    root = post.parent.root if (post and post.parent) else None

    # Initialize the renderer
    renderer = BiostarRenderer(escape=escape)

    # Initialize the lexer
    inline = BiostarInlineLexer(renderer=renderer, root=root, allow_rewrite=allow_rewrite)

    markdown = mistune.Markdown(hard_wrap=True, renderer=renderer, inline=inline)

    output = markdown(text=text)
    # Bleach clean the html.
    if clean:
        output = bleach.clean(text=output,
                            tags=ALLOWED_TAGS,
                            styles=ALLOWED_STYLES,
                            attributes=ALLOWED_ATTRIBUTES,
                            protocols=ALLOWED_PROTOCOLS)
    # Embed sensitive links into html
    output = linkify(text=output)

    return output


def test():
    html = parse(TEST_INPUT2)
    return html


if __name__ == '__main__':
    html = test()
    print(html)
