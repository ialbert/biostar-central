# Action codes.
BUMP_POST, OPEN_POST, TOGGLE_ACCEPT, MOVE_ANSWER, DUPLICATE, OFFTOPIC, DELETE, LOCK, CLOSE = range(9)

# Valid values for the order GET parameter.
RANK, VIEWS, REPLIES = ("rank", "views", "replies")

# The name of the session count data.
COUNT_DATA_KEY = "COUNT_DATA"

# Redirection field name.
REDIRECT_FIELD_NAME = 'next'

MYVOTES, MYPOSTS, MYTAGS, OPEN, FOLLOWING = ["myvotes", "myposts", "mytags", "open", "following"]

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
