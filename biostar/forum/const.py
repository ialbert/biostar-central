# Moderation action codes.
MOD_ACTIONS = list(range(9))
BUMP_POST, OPEN_POST, TOGGLE_ACCEPT, MOVE_ANSWER, DUPLICATE, OFFTOPIC, DELETE, CLOSE, REPORT_SPAM = MOD_ACTIONS

# Filter actions
RANK, VIEWS, REPLIES, TAGGED, \
VOTES, VISIT, REPUTATION, JOINED, \
ACTIVITY, UPDATE, ANSWERS, BOOKED, CREATION = ("rank", "views", "replies",
                    "tagged", "votes", "visit", "reputation",
                    "joined", "activity", "update", "answers", "bookmark", "creation")

# Map filter actions to respective database filters.
ORDER_MAPPER = {
    RANK: "-rank",
    UPDATE: "-lastedit_date",
    ANSWERS: "-answer_count",
    CREATION:"-creation_date",
    TAGGED: '-tagged',
    BOOKED: "-book_count",
    VIEWS: "-view_count",
    REPLIES: "-reply_count",
    VOTES: "-thread_votecount",
    VISIT: '-profile__last_login',
    REPUTATION: '-profile__score',
    JOINED: '-profile__date_joined',
    ACTIVITY: '-profile__date_joined'
}

LATEST_CACHE_KEY = "LATEST"

TAGS_CACHE_KEY = "TAGS"

SIMILAR_CACHE_KEY = "SIMILAR"
USERS_CACHE_KEY = "MENTIONED_USERS"

TAGS_FILE_KEY = "TAGS-FILE"

USERS_LIST_KEY = "USERS_LIST"
# The name of the session count data.
COUNT_DATA_KEY = "COUNT_DATA"

# Redirection field name.
REDIRECT_FIELD_NAME = 'next'

MYVOTES, MYPOSTS, MYTAGS, OPEN, FOLLOWING, SHOW_SPAM = ["myvotes", "myposts", "mytags", "open", "following", "spam"]

# Post type names for filtering in the UI.
JOBS, TOOLS, TUTORIALS, FORUM, PLANET = ["jobs", "tools", "tutorials", "forum", "planet"]

BOOKMARKS, MESSAGE = ["bookmarks", "message"]

MENTIONED, COMMUNITY, LATEST = ["mentioned", "community", "latest"]
