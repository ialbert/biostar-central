from biostar.settings import *


ALLOWED_TAGS = "p div br code pre h1 h2 h3 h4 hr span s sub sup b i img strong strike em underline super table "
ALLOWED_TAGS += "thead tr th td tbody"
ALLOWED_TAGS = ALLOWED_TAGS.split()

BUMP_POST, MOD_OPEN, TOGGLE_ACCEPT, MOVE_TO_ANSWER, \
MOVE_TO_COMMENT, DUPLICATE, CROSSPOST, CLOSE_OFFTOPIC, DELETE = range(9)

ALLOWED_STYLES = 'color font-weight background-color width height'.split()
ALLOWED_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],

}

REDIRECT_FIELD_NAME = 'next'

MYPOSTS, MYTAGS, OPEN, FOLLOWING = ["myposts", "mytags", "open", "following"]
BOOKMARKS, VOTES, MESSAGE = ["bookmarks", "votes", "message"]
MENTIONED, COMMUNITY, LATEST = ["mentioned", "community", "latest"]

# Topics that have a default tabs
TOPICS_WITH_TABS = [
                MYPOSTS, MYTAGS, FOLLOWING,
                BOOKMARKS, VOTES, MESSAGE, COMMUNITY, LATEST
              ]

PRIVATE_TOPICS = [ MYPOSTS, MYTAGS, MESSAGE, MENTIONED,
                    BOOKMARKS, FOLLOWING, VOTES
                    ]