# Action codes.
BUMP_POST, MOD_OPEN, TOGGLE_ACCEPT, MOVE_TO_ANSWER, \
MOVE_TO_COMMENT, DUPLICATE, CROSSPOST, CLOSE_OFFTOPIC, DELETE = range(9)

# Valid values for the order GET parameter.
RANK, VIEWS, REPLIES = ("rank", "views", "replies")

# Default tags shown on dropdowns
TAGS = ['RNA-Seq', 'ChIP-Seq', 'SNP', 'Assembly', 'software error', 'sequence',
             'subjective', 'general', 'solid', 'galaxy', 'motif', 'bed', 'conversion']


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
