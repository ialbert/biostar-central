"""
Various constants used throught the site

WARNING: DO NOT CHANGE ON PRODUCTION SEVERS! 
"""
import string

# The name of the editor group
MODERATOR_GROUP = 'mod_group'

# the minimal reputation needed to 
MIN_REP = 1

# post/user score change during an upvote
POST_SCORE_CHANGE = 1
USER_SCORE_CHANGE = 1

MAX_VOTES_PER_SESSION = 10    # this is how many votes can be cast per session

from datetime import timedelta
VOTE_SESSION_LENGTH = 60   # in seconds, the time intervals to reset vote limit
VOTE_SESSION_LENGTH = timedelta(seconds=VOTE_SESSION_LENGTH)

FIRST_SESSION = 'first-session'
LASTSORT_SESSION = 'last-sort'

# Add at the end
POST_QUESTION, POST_ANSWER, POST_COMMENT, POST_TUTORIAL, POST_BLOG, POST_FORUM, POST_NEWS, POST_REVIEW, POST_TOOL, POST_FIXME, POST_VIDEO, POST_JOB, POST_PUBLICATION, POST_TIP, POST_OTHER = range(1, 16)
POST_TYPES  = ( (POST_ANSWER, 'Answer') , (POST_COMMENT, 'Comment'), (POST_QUESTION, 'Question'), (POST_TUTORIAL, 'Tutorial'), (POST_TIP, 'Tip'),
    (POST_BLOG, 'Blog'), (POST_FORUM, 'Forum'), (POST_NEWS, 'News'), (POST_REVIEW, 'Review'), (POST_TOOL, 'Tool'), (POST_VIDEO, 'Video'),
    (POST_JOB, 'Job'), (POST_PUBLICATION, 'Research Paper') )

# direct mapping for quick lookups
POST_MAP  = dict( POST_TYPES )

# reverse mapping for quick lookups
POST_REV_MAP = dict( (y.lower(),x) for (x,y) in POST_MAP.items() )

# entities that will be displayed on the navigation bar
POST_NAV_BAR = [  ]
POST_NAV_BAR_LOWER = map(string.lower, POST_NAV_BAR)

# the valid sort orders
SORT_MAP = dict(
    rank="-rank", views="-views", creation="-creation_date",
    activity="-lastedit_date", votes="-full_score", answers="-answer_count",
    bookmark="-book_count", votes__date="-votes__date",
)

# valid pill entries
VALID_PILLS = set( "mytags all questions unanswered galaxy best bookmarked howto myvotes mybookmarks myposts jobs planet forum".split() )

# valid tab entries
VALID_TABS = set( "recent planet sticky".split() ) | VALID_PILLS 

# posts that only have content, no title or tags
POST_CONTENT_ONLY = set( [POST_ANSWER, POST_COMMENT ])

# these posts must have parent
POST_SUBLEVEL = set( [POST_ANSWER, POST_COMMENT ])

# main level posts 
POST_EXCLUDE = set( [POST_ANSWER, POST_COMMENT, POST_BLOG ])

# toplevel posts may stand alone and must have title and tags
POST_TOPLEVEL = set( POST_MAP.keys() ) - POST_SUBLEVEL

# posts the will go under forum
POST_FORUMLEVEL = set( (POST_FORUM, POST_NEWS, POST_REVIEW) )

# the session key that stores new post counts
SESSION_POST_COUNT = 'session-post-count'
SESSION_VIEW_COUNT = 'view-count'

# the type of messages that the system maintains
NOTE_USER, NOTE_MODERATOR, NOTE_ADMIN, NOTE_AWARD, NOTE_SITE, NOTE_PRIVATE = range(1, 7)
NOTE_TYPES = ((NOTE_USER,'User'), (NOTE_MODERATOR,'Moderator'), (NOTE_ADMIN,'Admin'), (NOTE_AWARD, 'Award'), (NOTE_SITE, "Site"), (NOTE_PRIVATE, "Private"))

# user types
USER_NEW, USER_MEMBER,  USER_MODERATOR, USER_ADMIN, USER_BLOG, USER_SPECIAL, USER_EXTERNAL = range(1, 8)
USER_TYPES = ( (USER_NEW, 'New'), (USER_MEMBER, 'Member'),  (USER_MODERATOR, 'Moderator'), 
    (USER_ADMIN, 'Administrator'), (USER_BLOG, 'Blog'), (USER_SPECIAL, 'Special'), (USER_EXTERNAL, "external"))

# user status types
USER_ACTIVE, USER_SUSPENDED, USER_BANNED = 1, 2, 3
USER_STATUS_TYPES = ( (USER_ACTIVE, 'Active'), (USER_SUSPENDED, 'Suspended'), (USER_BANNED, 'Banned' ))

# post status types        
POST_OPEN, POST_CLOSED, POST_DELETED = 100, 200, 300
POST_STATUS_TYPES = ( (POST_OPEN, 'Open'), (POST_CLOSED, 'Closed'), (POST_DELETED, 'Deleted') )

# the time between registering two post views
# from the same IP, in minutes
POST_VIEW_UPDATE = 30

# revision constants
REV_NONE, REV_CLOSE, REV_REOPEN, REV_DELETE, REV_UNDELETE = range(1000, 1005)
REV_ACTIONS = (
    (REV_NONE, ''), (REV_CLOSE, 'Close'), (REV_REOPEN, 'Reopen'),
    (REV_DELETE, 'Delete'), (REV_UNDELETE, 'Undelete')
)
REV_ACTION_MAP = dict(REV_ACTIONS)

# this stores the counts in the cache
CACHE_COUNT_KEY = "cache-count-key"

# moderation actions
USER_MODERATION, POST_MODERATION = 0, 1
USER_MOD_TYPES = [ (USER_MODERATION, 'Usermod'), (POST_MODERATION, 'Postmod') ]
    
# voting related constants
VOTE_UP, VOTE_DOWN, VOTE_ACCEPT, VOTE_BOOKMARK = range(1, 5)
VOTE_TYPES = ((VOTE_UP, 'Upvote'), (VOTE_DOWN, 'Downvote'), (VOTE_ACCEPT, 'Accept'), (VOTE_BOOKMARK, 'Bookmark'))
OPPOSING_VOTES = { VOTE_UP:VOTE_DOWN, VOTE_DOWN:VOTE_UP }

BADGE_BRONZE, BADGE_SILVER, BADGE_GOLD = 0, 1, 2
BADGE_TYPES = ((BADGE_BRONZE, 'bronze'), (BADGE_SILVER, 'silver'), (BADGE_GOLD, 'gold'))

BETA_TESTER_BADGE = "Beta Tester"

TARGET_COUNT_MAP = {
    POST_NEWS : "News",
    POST_QUESTION : "Question",
    POST_TOOL : "Tool",
    POST_TUTORIAL : "Tutorial",
    POST_JOB : "Job",
    POST_BLOG: "Blog",
    POST_VIDEO: "Video",
    "unanswered" : "Unanswered",
    POST_FORUM: "Forum",
}
