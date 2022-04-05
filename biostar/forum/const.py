
# Filter actions
VOTE_DATE, RANK, VIEWS, REPLIES, TAGGED, \
VOTES, VISIT, REPUTATION, JOINED, \
ACTIVITY, UPDATE, ANSWERS, BOOKED, CREATION = ('vote_date', "rank", "views", "replies",
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
    ACTIVITY: '-profile__date_joined',
    VOTE_DATE: '-votes__date'

}

ALLOWED_PARAMS = {"page", "order", "type", "limit", "query", "user", "active"}

# Cache keys used to cache objects.
LATEST_CACHE_KEY = "LATEST"
TAGS_CACHE_KEY = "TAGS"
SIMILAR_CACHE_KEY = "similar"
USERS_LIST_KEY = "USERS_LIST"

# The name of the session count data.
COUNT_DATA_KEY = "COUNT_DATA"
VOTES_COUNT = 'vote_count'

# Tabs to pick from in post listing
MYVOTES, MYPOSTS, MYTAGS, OPEN, \
FOLLOWING, SHOW_SPAM, BOOKMARKS = ["myvotes", "myposts", "mytags", "open", "following", "spam", "bookmarks"]

JOBS, TOOLS, TUTORIALS, FORUM, PLANET, LATEST = ["jobs", "tools", "tutorials", "forum", "planet", "latest"]

