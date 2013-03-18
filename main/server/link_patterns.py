# coding=utf-8
"""
Link patterns matched in the markdown2 library.

Requires a modified markdown2 library!

"""

__author__ = 'ialbert'

from django.conf import settings
import re

DOMAIN = settings.SITE_DOMAIN

class AutoLink(object):
    pattern = re.compile(r'((?P<start>(\s|\(|^))(?P<url>((http://|ftp://|https://)\S+)))', re.I)

    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
        url  = m.group('url')
        start = m.group('start')
        link = start + "<a href='%s'>%s</a>"

        # some common corner cases matched explicitly rather than regexp matching to simplify the RE
        # recognize various punctuations at the end
        if url[-1] in r":,.()'\<>\"":
            end = url[-1]
            url = url[:-1]
            link = link + end

        return link % (url, url)

class TagMatch(object):

    pattern1 = re.compile(r"(?P<url>(http://%s/show/tag/)(?P<tag>(\S+)))/" % DOMAIN, re.I)
    pattern2 = re.compile(r"(tag\s*?:\s*?(?P<tag>(\S+)))", re.I)

    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
        tag = m.group('tag')
        link = "tag: <a href='/show/tag/%s/'>%s</a>" % (tag, tag)
        return link

class UserMatch(object):

    pattern1 = re.compile(r"(?P<url>(http://%s/u/)(?P<uid>(\d+)))/" % DOMAIN, re.I)
    pattern2 = re.compile(r"(?P<url>(\\user\s)(?P<uid>(\d+)))", re.I)
    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
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

class PostMatch(object):

    pattern1 = re.compile(r"(?P<url>(http://%s/p/)(?P<uid>(\S+)))/" % DOMAIN, re.I)
    pattern2 = re.compile(r"(?P<url>(\\post\s+)(?P<uid>(\d+)))", re.I)
    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
        from main.server import models
        uid = m.group('uid')
        url = m.group('url')
        try:
            posts = models.Post.objects.filter(id=uid)
            post = posts[0]
            link = "<a href='%s'>%s</a>" % (post.get_absolute_url(), post.title)
        except Exception, exc:
            print "*** postmatch error %s, %s, %s" %(url, uid, exc)
            link = "%s%s" % (url, uid)

        return link

class GistMatch(object):
    pattern1 = re.compile(r"(\\gist\s+(?P<uid>(\d+)))", re.I)
    pattern2 = re.compile(r"(gist\s*?:\s*?(?P<uid>(\d+)))", re.I)
    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
        uid = m.group('uid')
        link = "<div><script src='http://gist.github.com/%s.js'></script></div>" % uid
        return link

class YouTubeMatch(object):
    pattern1 = re.compile(r"(\\youtube\s+?(?P<uid>(\S+)))", re.I)
    pattern2 = re.compile(r"(youtube\s*?:\s*?(?P<uid>(\S+)))", re.I)
    @staticmethod
    def action(m):
        "This matches a user based on a user pattern"
        uid = m.group('uid')
        link = r'''
        <div>
            <iframe width='560' height='315' src='http://www.youtube.com/embed/%s' frameborder='0' allowfullscreen></iframe>
        </div>
        ''' % uid
        return link

classes = [
    TagMatch, YouTubeMatch, GistMatch, PostMatch, UserMatch
]
all = []
for c in classes:
    all.append((c.pattern1, c.action))
    all.append((c.pattern2, c.action))

# add the autolink
all.append((AutoLink.pattern, AutoLink.action))

