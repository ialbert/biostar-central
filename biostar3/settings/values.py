#
# Site specic behaviors.
#
from collections import OrderedDict

# How many recent votes to show.
RECENT_VOTE_COUNT = 10

# How many recent users to show.
RECENT_USER_COUNT = 10

# How many posts per page.
POST_PAGINATE_BY = 10

# The CSS classes associated with the Django messages framework.
MESSAGE_TAGS = {
    10: 'alert-info', 20: 'alert-info', 25: 'alert-success', 30: 'alert-warning', 40: 'alert-danger',
}

# Connects a word in post sort to a queryset sort attribute of the data model.
POST_SORT_MAP = OrderedDict([
    ("update", "-lastedit_date"),
    ("views", "-view_count"),
    ("followers", "-subs_count"),
    ("answers", "-reply_count"),
    ("bookmarks", "-book_count"),
    ("votes", "-vote_count"),
    ("rank", "-rank"),
    ("creation", "-creation_date"),
])

# These are the fields rendered in the post sort order drop down.
POST_SORT_FIELDS = POST_SORT_MAP.keys()
POST_SORT_DEFAULT = POST_SORT_FIELDS[0]

# The messages show when the sort is not valid.
POST_SORT_INVALID_MSG = "Invalid sort parameter received"