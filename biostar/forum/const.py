# Action codes.
BUMP_POST, OPEN_POST, TOGGLE_ACCEPT, MOVE_ANSWER, DUPLICATE, OFFTOPIC, DELETE, REPORT_SPAM = range(8)

# Valid values for the order GET parameter.
RANK, VIEWS, REPLIES = ("rank", "views", "replies")

MYVOTES_CACHE_KEY = "MYVOTES"
TAGS_CACHE_KEY = "TAGS"
MYPOSTS_CACHE_KEY = "MYPOSTS"
FOLLOWING_CACHE_KEY = "FOLLOWING"
BOOKMARKS_CACHE_KEY = "BOOKMARKS"
MYTAGS_CACHE_KEY = "MYTAGS"

SIMILAR_CACHE_KEY = "SIMILAR"
USERS_CACHE_KEY = "MENTIONED_USERS"


# The name of the session count data.
COUNT_DATA_KEY = "COUNT_DATA"

# Redirection field name.
REDIRECT_FIELD_NAME = 'next'

MYVOTES, MYPOSTS, MYTAGS, OPEN, FOLLOWING, SHOW_SPAM = ["myvotes", "myposts", "mytags", "open", "following", "spam"]

BOOKMARKS, MESSAGE = ["bookmarks", "message"]

MENTIONED, COMMUNITY, LATEST = ["mentioned", "community", "latest"]

# Post type names for filtering in the UI.
JOBS, TOOLS, TUTORIALS, FORUM, PLANET = ["jobs", "tools", "tutorials", "forum", "planet"]

# Topics that have a default tabs
TOPICS_WITH_TABS = [
    MYPOSTS, MYTAGS, FOLLOWING,
    BOOKMARKS, MYVOTES, MESSAGE, COMMUNITY, LATEST
]

PRIVATE_TOPICS = [MYPOSTS, MYTAGS, MESSAGE, MENTIONED,
                  BOOKMARKS, FOLLOWING, MYVOTES
                  ]
