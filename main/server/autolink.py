# coding=utf-8
"""

Requires a modified markdown2 library!

"""
__author__ = 'ialbert'

from django.conf import settings
from django.core.urlresolvers import reverse

import re

DOMAIN = settings.SITE_DOMAIN

BIOINFO_WORDS = re.compile(r"(?P<word>\b(bwa|sam|bam|samtools|bedtools|sam)\b)", re.I)

BIOINFO_PATT = {
    'bwa' : 'http://bio-bwa.sourceforge.net/',
    'sam' : 'http://samtools.sourceforge.net/SAM1.pdf',
    'bam' : 'http://samtools.sourceforge.net/SAM1.pdf',
    'samtools' : 'http://samtools.sourceforge.net/',
    'bedtools' : 'https://code.google.com/p/bedtools/',
}
def bioinfo_link(m):
    word = m.group('word')
    patt = BIOINFO_PATT.get(word.lower(),"?")
    link = "<a href='%s'>%s</a>" % (patt, word)
    return link

TAG_FULL = re.compile(r"(?P<url>http://%s/show/tag/)(?P<tag>\S+)/" % DOMAIN, re.I)

def tag_link(m):
    "Tag link generator"
    tag = m.group('tag')
    link = "<a href='%s'>%s</a>" % (reverse("show-tag", args=[tag]), tag)
    return link

USER_FULL = re.compile(r"(?P<url>http://%s/u/)(?P<uid>\d+)/*" % DOMAIN, re.I)
USER_SHORT = re.compile(r"(?P<url>\\user\s)(?P<uid>\d+)", re.I)

def user_link(m):
    "User link generator"
    from main.server import models

    uid = m.group('uid')
    url = m.group('url')
    users = models.User.objects.filter(id=uid).select_related('user__profile')
    if users:
        user = users[0]
        link = "<a href='%s'>%s</a>" % (user.profile.get_absolute_url(), user.profile.display_name)
    else:
        link = "%s%s" % (url, uid)

    return link

POST_FULL = re.compile(r"(?P<url>http://%s/p/\d+/#)(?P<uid>\d+)" % DOMAIN, re.I)
POST_TOP  = re.compile(r"(?P<url>http://%s/p/)(?P<uid>\d+)/*" % DOMAIN, re.I)
POST_SHORT = re.compile(r"(?P<url>\\post\s+)(?P<uid>\d+)", re.I)

def post_link(m):
    "Post link generator"
    from main.server import models

    uid = m.group('uid')
    url = m.group('url')
    posts = models.Post.objects.filter(id=uid)
    if posts:
        post = posts[0]
        link = "<a href='%s'>%s</a>" % (post.get_absolute_url(), post.title)
    else:
        link = "%s%s" % (url, uid)

    return link


GIST_FULL = re.compile(r"(?P<url>https://gist.github.com/)(?P<uid>\S+)", re.I)
GIST_SHORT = re.compile(r"(?P<url>(\\gist\s+)(?P<uid>(\S+)))", re.I)


def gist_link(m):
    "Creates Gist links"

    uid = m.group('uid')
    url = m.group('url')

    if url.startswith("\\"):
        url = "https://gist.github.com/"

    # trim off extra js
    if uid.endswith(".js"):
        uid = uid[:-4]

    link = "<div><script src='%s%s.js'></script></div>" % (url, uid)
    return link


YOUTUBE_FULL = re.compile(r"(?P<url>http://www.youtube.com/watch\?v=)(?P<uid>\S+)", re.I)
YOUTUBE_SHORT = re.compile(r"(?P<url>\\youtube\s+)(?P<uid>\S+)", re.I)


def youtube_link(m):
    "Creates YouTube links"
    url = m.group("url")
    uid = m.group('uid')
    link = r'''
        <div>
            <iframe width='560' height='315' src='http://www.youtube.com/embed/%s' frameborder='0' allowfullscreen></iframe>
        </div>
        ''' % uid
    return link


AUTO_LINK_PATTERN = re.compile(r'(\s|^)(?P<url>(http|https|ftp)://\S+)', re.I)
EXCLUDE_CHARS = ";,)].:"


def auto_link(m):
    "Auto links"
    url = m.group('url')
    end = ""

    # a "real" URL regexp is very complex so we do it this way
    # these are characters that are usually a punctuation at the end
    while 1:
        if url[-1] in EXCLUDE_CHARS:
            end += url[-1]
            url = url[:-1]
            continue
        break
    end = end[::-1]
    link = "<a href='%s'>%s</a>%s" % (url, url, end)
    return link

patterns = [

    (BIOINFO_WORDS, bioinfo_link),

    (TAG_FULL, tag_link),

    (USER_SHORT, user_link),
    (USER_FULL, user_link),

    (POST_SHORT, post_link),
    (POST_FULL, post_link),
    (POST_TOP, post_link),

    (GIST_FULL, gist_link),
    (GIST_SHORT, gist_link),

    (YOUTUBE_FULL, youtube_link),
    (YOUTUBE_SHORT, youtube_link),

    (AUTO_LINK_PATTERN, auto_link),

]

def test():
    text = """

(http://www.psu.edu;

 - (http://www.psu.edu).

(http://www.psu.edu) http://www.psu.edu. ok


\gist ialbert/ae46c5f51d63cdf2d0d2

https://gist.github.com/ialbert/ae46c5f51d63cdf2d0d2

\post 98

http://%(domain)s/p/22/#98

\user 2

http://%(domain)s/u/2/

http://%(domain)s/show/tag/bwa+galaxy+bowtie/

Youtube Embedding
-----------------

\youtube tY2n2CHMXfI

http://www.youtube.com/watch?v=tY2n2CHMXfI


    """ % dict(domain=DOMAIN)

    import html

    page = html.generate(text)

    print page


if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    test()



