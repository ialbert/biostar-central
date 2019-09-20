"""
Markdown parser to render the Biostar style markdown.
"""
import re

import mistune
import requests
from django.shortcuts import reverse
from django.db.models import F
from django.conf import settings
from mistune import Renderer, InlineLexer, InlineGrammar

from biostar.forum import auth
from biostar.forum.models import Post, Subscription
from biostar.accounts.models import Profile, User

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
SITE_URL = f"{settings.SITE_DOMAIN}{settings.HTTP_PORT}"


# Biostar patterns
PORT = ':' + settings.HTTP_PORT if settings.HTTP_PORT else ''

USER_PATTERN = rec(fr"^http(s)?://{settings.SITE_DOMAIN}{PORT}/accounts/profile/(?P<uid>[\w_.-]+)(/)?")
POST_TOPLEVEL = rec(fr"^http(s)?://{settings.SITE_DOMAIN}{PORT}/p/(?P<uid>(\w+))(/)?$")
POST_ANCHOR = rec(fr"^http(s)?://{settings.SITE_DOMAIN}{PORT}/p/\w+/\#(?P<uid>(\w+))(/)?")
# Match any alphanumeric characters after the @.
# These characters are allowed in handles: _  .  -
MENTINONED_USERS = rec(r"(\@(?P<handle>[\w_.-]+))")

# Youtube pattern.
YOUTUBE_PATTERN1 = rec(r"^http(s)?://www.youtube.com/watch\?v=(?P<uid>([\w-]+))(/)?")
YOUTUBE_PATTERN2 = rec(r"https://www.youtube.com/embed/(?P<uid>([\w-]+))(/)?")
YOUTUBE_PATTERN3 = rec(r"https://youtu.be/(?P<uid>([\w-]+))(/)?")
YOUTUBE_HTML = '<iframe width="420" height="315" src="//www.youtube.com/embed/%s" frameborder="0" allowfullscreen></iframe>'

# Ftp link pattern.
FTP_PATTERN = rec(r"^ftp://[\w\.]+(/?)$")

# Gist pattern.
GIST_PATTERN = rec(r"^https://gist.github.com/(?P<uid>([\w/]+))")
GIST_HTML = '<script src="https://gist.github.com/%s.js"></script>'

# Twitter pattern.
TWITTER_PATTERN = rec(r"http(s)?://twitter.com/\w+/status(es)?/(?P<uid>([\d]+))")


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
    """
    text = re.compile(r'^[\s\S]+?(?=[\\<!\[*`~@]|https?://| {2,}\n|$)')


class BiostarInlineLexer(MonkeyPatch):
    grammar_class = BiostarInlineGrammer

    def __init__(self, root=None, *args, **kwargs):
        self.root = root

        super(BiostarInlineLexer, self).__init__(*args, **kwargs)

    def enable_post_link(self):
        self.rules.post_link = POST_TOPLEVEL
        self.default_rules.insert(0, 'post_link')

    def enable_mention_link(self):
        self.rules.mention_link = MENTINONED_USERS
        self.default_rules.insert(0, 'mention_link')

    def output_mention_link(self, m):

        handle = m.group("handle")
        # Query user and get the link
        user = User.objects.filter(username=handle).first()
        if user:
            profile = reverse("user_profile", kwargs=dict(uid=user.profile.uid))
            link = f'<a href="{profile}">{user.profile.name}</a>'
            # Subscribe mentioned users to post.
            if self.root:
                # Create user subscription if it does not already exist.
                auth.create_subscription(post=self.root, user=user, delete_exisiting=False)
        else:
            link = m.group(0)

        return link

    def output_post_link(self, m):
        uid = m.group("uid")
        post = Post.objects.filter(uid=uid).first() or Post(title=f"Invalid post uid: {uid}")
        link = m.group(0)
        return f'<a href="{link}">{post.title}</a>'

    def enable_anchor_link(self):
        self.rules.anchor_link = POST_ANCHOR
        self.default_rules.insert(0, 'anchor_link')

    def output_anchor_link(self, m):
        uid = m.group("uid")
        alt, link = f"{uid}", m.group(0)
        return f'<a href="{link}">ANCHOR: {alt}</a>'

    def enable_user_link(self):
        self.rules.user_link = USER_PATTERN
        self.default_rules.insert(0, 'user_link')

    def output_user_link(self, m):
        uid = m.group("uid")
        link = m.group(0)
        profile = Profile.objects.filter(uid=uid).first()
        name = profile.name if profile else f"Invalid user uid: {uid}"

        return f'<a href="{link}">USER: {name}</a>'

    def enable_youtube_link1(self):
        self.rules.youtube_link1 = YOUTUBE_PATTERN1
        self.default_rules.insert(1, 'youtube_link1')

    def output_youtube_link1(self, m):
        uid = m.group("uid")
        return YOUTUBE_HTML % uid

    def enable_youtube_link2(self):
        self.rules.youtube_link2 = YOUTUBE_PATTERN2
        self.default_rules.insert(1, 'youtube_link2')

    def output_youtube_link2(self, m):
        uid = m.group("uid")
        return YOUTUBE_HTML % uid

    def enable_youtube_link3(self):
        self.rules.youtube_link3 = YOUTUBE_PATTERN3
        self.default_rules.insert(1, 'youtube_link3')

    def output_youtube_link3(self, m):
        uid = m.group("uid")
        return YOUTUBE_HTML % uid

    def enable_twitter_link(self):
        self.rules.twitter_link = TWITTER_PATTERN
        self.default_rules.insert(1, 'twitter_link')

    def output_twitter_link(self, m):
        uid = m.group("uid")
        return get_tweet(uid)

    def enable_gist_link(self):
        self.rules.gist_link = GIST_PATTERN
        self.default_rules.insert(3, 'gist_link')

    def output_gist_link(self, m):
        uid = m.group("uid")
        return GIST_HTML % uid

    def enable_ftp_link(self):
        self.rules.ftp_link = FTP_PATTERN
        self.default_rules.insert(3, 'ftp_link')

    def output_ftp_link(self, m):
        link = m.group(0)
        return f'<a href="{link}">{link}</a>'


def parse(text, post=None):
    """
    Parses markdown into html.
    Expands certain patterns into HTML.
    """
    # Resolve the root
    root = post.parent.root if (post and post.parent) else None
    inline = BiostarInlineLexer(renderer=Renderer(), root=root)
    inline.enable_post_link()
    inline.enable_mention_link()
    inline.enable_anchor_link()
    inline.enable_anchor_link()
    inline.enable_user_link()

    inline.enable_youtube_link1()
    inline.enable_youtube_link2()
    inline.enable_youtube_link3()
    inline.enable_ftp_link()
    inline.enable_twitter_link()

    markdown = mistune.Markdown(escape=True, hard_wrap=True, inline=inline)

    html = markdown(text)

    return html


def test():
    html = parse(TEST_INPUT2)
    return html


if __name__ == '__main__':

    html = test()
    print(html)
