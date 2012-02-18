"""
Various constants used throught the site

WARNING: DO NOT CHANGE ON PRODUCTION SEVERS! 
"""
# The name of the editor group
MODERATOR_GROUP = 'mod_group'

# the minimal reputation needed to 
MIN_REP = 1

# post/user score change during an upvote
POST_SCORE_CHANGE = 1
USER_SCORE_CHANGE = 1

VOTE_SESSION_LENGTH   = 60  # in seconds
MAX_VOTES_PER_SESSION = 3   # this is how many votes can be cast per session

# Add at the end
POST_QUESTION, POST_ANSWER, POST_COMMENT, POST_GUIDE, POST_BLOG, POST_NEWS, POST_FORUM, POST_RANT, POST_REVIEW, POST_OTHER = range(1, 11)
POST_TYPES  = ( (POST_ANSWER, 'Answer') , (POST_COMMENT, 'Comment'), (POST_QUESTION, 'Question'), (POST_GUIDE, 'Guide'), 
    (POST_NEWS, 'News'), (POST_BLOG, 'Blog'), (POST_FORUM, 'Forum'), (POST_REVIEW, 'Review'), (POST_RANT, 'Rant'))

# direct mapping for quick lookups
POST_MAP  = dict( POST_TYPES )

# reverse mapping for quick lookups
POST_REV_MAP = dict( (y.lower(),x) for (x,y) in POST_MAP.items() )

# posts that only have content, no title or tags
POST_CONTENT_ONLY = set( [POST_ANSWER, POST_COMMENT ])

# these posts must have parent
POST_SUBLEVEL = set( [POST_ANSWER, POST_COMMENT ])

# toplevel posts may stand alone and must have title and tags
POST_TOPLEVEL = set( POST_MAP.keys() ) - POST_SUBLEVEL

# the session key that stores new post counts
SESSION_POST_COUNT = 'session-post-count'

# the type of messages that the system maintains
NOTE_USER, NOTE_MODERATOR, NOTE_ADMIN, NOTE_AWARD, NOTE_SITE = range(1, 6)
NOTE_TYPES = ((NOTE_USER,'User'), (NOTE_MODERATOR,'Moderator'), (NOTE_ADMIN,'Admin'), (NOTE_AWARD, 'Award'), (NOTE_SITE, "Site"))

# user types
USER_NORMAL,  USER_MODERATOR, USER_ADMIN, USER_SPECIAL, = range(1, 5)
USER_TYPES = ( (USER_NORMAL, 'Member'),  (USER_MODERATOR, 'Moderator'), 
    (USER_ADMIN, 'Administrator'), (USER_SPECIAL, 'Special'),)

# user status types
USER_ACTIVE, USER_SUSPENDED = 10, 20
USER_STATUS_TYPES = ( (USER_ACTIVE, 'Active'), (USER_SUSPENDED, 'Suspended') )

# post status types        
POST_OPEN, POST_CLOSED, POST_DELETED = 100, 200, 300
POST_STATUS_TYPES = ( (POST_OPEN, 'Open'), (POST_CLOSED, 'Closed'), (POST_DELETED, 'Deleted') )

# revision constants
REV_NONE, REV_CLOSE, REV_REOPEN, REV_DELETE, REV_UNDELETE = range(1000, 1005)
REV_ACTIONS = (
    (REV_NONE, ''), (REV_CLOSE, 'Close'), (REV_REOPEN, 'Reopen'),
    (REV_DELETE, 'Delete'), (REV_UNDELETE, 'Undelete')
)
REV_ACTION_MAP = dict(REV_ACTIONS)

# moderation actions
USER_MODERATION, POST_MODERATION = 0, 1
USER_MOD_TYPES = [ (USER_MODERATION, 'Usermod'), (POST_MODERATION, 'Postmod') ]
    
# voting related constants
VOTE_UP, VOTE_DOWN, VOTE_ACCEPT, VOTE_BOOKMARK = range(1, 5)
VOTE_TYPES = ((VOTE_UP, 'Upvote'), (VOTE_DOWN, 'Downvote'), (VOTE_ACCEPT, 'Accept'), (VOTE_BOOKMARK, 'Bookmark'))

BADGE_BRONZE, BADGE_SILVER, BADGE_GOLD = 0, 1, 2
BADGE_TYPES = ((BADGE_BRONZE, 'bronze'), (BADGE_SILVER, 'silver'), (BADGE_GOLD, 'gold'))

BETA_TESTER_BADGE = "Beta Tester"
