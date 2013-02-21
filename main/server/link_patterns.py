"""
Link patterns to process the markdown
"""

__author__ = 'ialbert'

import re

def code_match(m):
    value = m.group('value')
    return r"http://code.activestate.com/recipes/%s/" % value

code_match.patt = re.compile("recipe\s+(?P<id>(\d+))", re.I)


class UserMatch(object):
    pattern = re.compile("(?P<url>(http://\S+?/u/)(?P<uid>(\d+)))/", re.I)

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


all = [
    (code_match.patt, code_match),
    (UserMatch.pattern, UserMatch.action),
    ]