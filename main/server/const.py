"""
Various constants used throught the site

WARNING: DO NOT CHANGE ON PRODUCTION SEVERS! 
"""

# The name of the editor group
MODERATOR_GROUP = 'mod_group'

# the minimal reputation needed to 
MIN_REP = 1

# Permissions granted only to admins and moderators. Admins have mod permissions.
MODERATOR_PERM = ['moderate_post', 'view_deleted']
ADMIN_PERM = []

# Add at the end
POST_QUESTION, POST_ANSWER, POST_COMMENT, POST_GUIDE, POST_BLOG, POST_NEWS, POST_OTHER = range(1, 8)
POST_TYPES  = ( (POST_QUESTION, 'Question'), (POST_GUIDE, 'Guide'), (POST_BLOG, 'Blog'), 
    (POST_NEWS, 'News'), (POST_OTHER, 'Post'), (POST_ANSWER, 'Answer') , (POST_COMMENT, 'Comment'),  )

# for quick lookups
POST_MAP  = dict( POST_TYPES )

# posts that require full form
POST_FULL_FORM = set( (POST_QUESTION, POST_GUIDE, POST_BLOG) )

# the type of messages that the system maintains
NOTE_USER, NOTE_MODERATOR, NOTE_ADMIN, NOTE_AWARD, NOTE_SITE = range(1, 6)
NOTE_TYPES = ((NOTE_USER,'User'), (NOTE_MODERATOR,'Moderator'), (NOTE_ADMIN,'Admin'), (NOTE_AWARD, 'Award'), (NOTE_SITE, "Site"))

# user types
USER_NORMAL,  USER_MODERATOR, USER_ADMIN, USER_SPECIAL, = range(1, 5)
USER_TYPES = ( (USER_NORMAL, 'Member'),  (USER_MODERATOR, 'Moderator'), 
    (USER_ADMIN, 'Administrator'), (USER_SPECIAL, 'Special'),)

# revision constants
REV_NONE, REV_CLOSE, REV_REOPEN, REV_DELETE, REV_UNDELETE = range(1, 6)
REV_ACTIONS = (
    (REV_NONE, ''), (REV_CLOSE, 'Close'), (REV_REOPEN, 'Reopen'),
    (REV_DELETE, 'Delete'), (REV_UNDELETE, 'Undelete')
)
REV_ACTION_MAP = dict(REV_ACTIONS)

# moderation actions
USER_MODERATION, POST_MODERATION = 0, 1
USER_MOD_CHOICES = [ (USER_MODERATION, 'Usermod'), (POST_MODERATION, 'Postmod') ]
    
# voting related constants
VOTE_UP, VOTE_DOWN, VOTE_ACCEPT, VOTE_FAVORITE = range(1, 5)
VOTE_TYPES = ((VOTE_UP, 'Upvote'), (VOTE_DOWN, 'Downvote'), (VOTE_ACCEPT, 'Accept'), (VOTE_FAVORITE, 'Favorite'))

# mappings of mutually exclusive votes
OPPOSING_VOTES = { VOTE_UP:VOTE_DOWN, VOTE_DOWN:VOTE_UP } 

# post score changes
POST_SCORE = { VOTE_UP:1, VOTE_DOWN:-1, VOTE_FAVORITE:2 }

# user reputation changes
USER_REP  = {
    VOTE_UP:10,
    VOTE_DOWN:-2,
    VOTE_ACCEPT:15
}

# voter reputation changes
VOTER_REP = {
    VOTE_DOWN: -1,
    VOTE_ACCEPT:2
}

BADGE_BRONZE, BADGE_SILVER, BADGE_GOLD = 0, 1, 2
BADGE_TYPES = ((BADGE_BRONZE, 'bronze'), (BADGE_SILVER, 'silver'), (BADGE_GOLD, 'gold'))

BETA_TESTER_BADGE = "Beta Tester"