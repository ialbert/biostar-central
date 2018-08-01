from biostar.settings import *


ALLOWED_TAGS = "p div br code pre h1 h2 h3 h4 hr span s sub sup b i img strong strike em underline super table "
ALLOWED_TAGS += "thead tr th td tbody"
ALLOWED_TAGS = ALLOWED_TAGS.split()

ALLOWED_STYLES = 'color font-weight background-color width height'.split()
ALLOWED_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],

}

REDIRECT_FIELD_NAME = 'next'
TOPIC_FIELD_NAME = "topic"
MESSAGE_TABS = ["inbox", "outbox", "unread"]

INBOX, OUTBOX, UNREAD = MESSAGE_TABS

ACTIVE_TAB = "active"

MYPOSTS, MYTAGS, UNANSWERED, FOLLOWING = ["myposts", "mytags", "open", "following"]
BOOKMARKS, VOTES, MESSAGE = ["bookmarks", "votes", "message"]
MENTIONED, COMMUNITY, LATEST = ["mentioned", "community", "latest"]


PRIVATE_TOPICS = [ MYPOSTS, MYTAGS, MESSAGE, MENTIONED, OUTBOX,
                    INBOX, BOOKMARKS, UNREAD, FOLLOWING, VOTES
                    ]