"""
Various constants used throught the site
"""

# The name of the editor group
MODERATOR_GROUP = 'mod_group'

# the minimal reputation needed to 
MIN_REP = 1


# Permissions granted only to admins and moderators. Admins have mod permissions.
MODERATOR_PERM = ['moderate_post', 'view_deleted']
ADMIN_PERM = []

# I added here anything that I could think of for the future rather than
# all post types that are to be implemented. This could make life a lot easier in the future
POST_NAMES =  "Post Question Answer Comment Guide Blog News Opinion Announcement Article Pony Other".split()
POST_VALS  = range(len(POST_NAMES))

# this can go into the models
POST_CHOICES  = zip( POST_VALS, POST_NAMES )

# allows quick check for valid post names or values
POST_MAP      = dict(zip( POST_NAMES, POST_VALS ))
POST_REV_MAP  = dict( [ (v, k) for k,v in POST_MAP.items() ] )

# convenience constants
POST_QUESTION = POST_MAP['Question']
POST_ANSWER   = POST_MAP['Answer']
POST_COMMENT  = POST_MAP['Comment']

# user types
USER_NORMAL, USER_MODERATOR, USER_ADMIN = range(0, 3)
USER_TYPES = ((USER_NORMAL, 'Member'), (USER_MODERATOR, 'Moderator'), (USER_ADMIN, 'Administrator'))

# revision constants
REV_NONE, REV_CLOSE, REV_REOPEN, REV_DELETE, REV_UNDELETE = range(0, 5)
REV_ACTIONS = (
    (REV_NONE, ''), (REV_CLOSE, 'Close'), (REV_REOPEN, 'Reopen'),
    (REV_DELETE, 'Delete'), (REV_UNDELETE, 'Undelete')
)
REV_ACTION_MAP = dict(REV_ACTIONS)

# moderation actions
USER_MODERATION, POST_MODERATION = 0, 1
USER_MOD_CHOICES = [ (USER_MODERATION, 'Usermod'), (POST_MODERATION, 'Postmod') ]
    
# voting related constants
VOTE_UP, VOTE_DOWN, VOTE_ACCEPT = range(0, 3)
VOTE_TYPES = ((VOTE_UP, 'Upvote'), (VOTE_DOWN, 'Downvote'), (VOTE_ACCEPT, 'Accept'))

# mappings of mutually exclusive votes
OPPOSING_VOTES = { VOTE_UP:VOTE_DOWN, VOTE_DOWN:VOTE_UP } 

# post score changes
POST_SCORE = { VOTE_UP:1, VOTE_DOWN:-1 }

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

