from biostar.settings import *


BUMP_POST, MOD_OPEN, TOGGLE_ACCEPT, MOVE_TO_ANSWER, \
MOVE_TO_COMMENT, DUPLICATE, CROSSPOST, CLOSE_OFFTOPIC, DELETE = range(9)

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