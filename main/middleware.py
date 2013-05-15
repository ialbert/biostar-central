from datetime import datetime, timedelta
from itertools import *
from django.contrib.auth import logout
from django.contrib import messages
from django.conf import settings
from main.server import models, notegen, html, awards
from main.server.const import *
from django.core.cache import cache

settings.CONTENT_INDEXING = True

class Session(object):
    """
    This object maintains various session data for a user. To minimize database access
    saves all at just one time. It also works for non-authenticated users but avoids
    creating a database sessions for them
    """
    SESSION_KEY, SORT_KEY, COUNT_KEY, TAB, ALL, RANK = "session-data", 'sortz', 'count', 'tabz', 'all', 'activity'
    def __init__(self, request):
        self.request = request
        self.has_storage = request.user.is_authenticated()
        default = { self.COUNT_KEY:{ }, self.SORT_KEY:self.RANK, self.TAB:self.ALL}
        if self.has_storage:
            self.data = self.request.session.get(self.SESSION_KEY, default )
        else:
            self.data = default
    
    def save(self):
        "Saves the counts back to the session"
        if self.has_storage:
            self.request.session[self.SESSION_KEY] = self.data
    
    def set_tab(self, value):
        self.data[self.TAB] = value
        
    def get_tab(self, value=None):
        "Facilitates navigation by remebering the last visited url"    
        return self.data.get(self.TAB, self.ALL)
    
    def sort_order(self):
        "Stores the last sort order in the session"
        value = self.request.GET.get('sort', '').lower()
        value = value or self.data.get(self.SORT_KEY)
        if value not in SORT_MAP:
            value = self.RANK
        # disable storing sort order in the sessions
        #self.data[self.SORT_KEY] = value
        return value

    def set_counts(self, value):
        self.data[self.COUNT_KEY] = value
            
    def get_counts(self, post_type=None):
        if not self.has_storage:
            return generate_counts(self.request)
        else:

            key = TARGET_COUNT_MAP.get(post_type, None)
            self.data[self.COUNT_KEY][key] = 0
            return self.data[self.COUNT_KEY]
            
def generate_counts(request, weeks=settings.COUNT_INTERVAL_WEEKS):
    "Returns the number of counts for each post type in the interval that has passed"
    user = request.user
    now  = datetime.now()

    counts = cache.get(CACHE_COUNT_KEY)

    if counts:
        return counts

    if user.is_authenticated():
        since = user.profile.last_visited
    else:
        since = now - timedelta(weeks=weeks)

    # for debugging
    #since = now - timedelta(weeks=1000)

    # posts since the last visit
    pairs = models.Post.objects.filter(type__in=POST_TOPLEVEL, status=POST_OPEN, creation_date__gt=since).order_by('-id').values_list("type", "answer_count")
    
    # establish how many of the posts have not been answered 
    values = [ p[0] for p in pairs ]
    unanswered = len([ ptype for (ptype, pcount) in pairs if (ptype == POST_QUESTION) and (pcount == 0) ])
    
    # needs to be sorted for the groupby
    values.sort()

    # how many times does each post type appear in the list
    counts = dict( [ (POST_MAP[k], len(list(v))) for (k, v) in groupby(values) ] )
    
    # fill in unanswered posts
    counts['Unanswered'] = unanswered

    counts['howto'] = counts.get("Tutorial", 0) + counts.get("Tool", 0) + counts.get("Tip", 0)

    if user.is_authenticated():
        vote_count = models.Vote.objects.filter(post__author=user, date__gt=since).count()
        counts['vote_count'] = vote_count
        counts['message_count'] = user.profile.new_messages

    if not user.is_authenticated():
        # store the cache key for non-authenticated users
        cache.set(CACHE_COUNT_KEY, counts, 3600)

    return counts

class LastVisit(object):
    """
    Updates the last visit stamp at SESSION_UPDATE_TIME intervals
    """

    def process_request(self, request):
        user = request.user
        sess = Session(request)
        
        # check suspended status for users
        if user.is_authenticated() and ( user.profile.suspended ):
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
            return
            
        # handle anonymous users
        if not user.is_authenticated():
            # anonymous users
            request.user.can_moderate = False
            return 
            
        # at this point we only have authenticated users
        profile = user.get_profile()
        
        # setting a handy shortcut
        request.user.can_moderate = profile.can_moderate
        
        # only write to database intermittently
        expired = (datetime.now() - profile.last_visited).seconds
            
        if expired > settings.SESSION_UPDATE_TIME:
            counts = generate_counts(request)
            sess.set_counts(counts)
            sess.save()

            # save the last update time
            profile.update_expiration()

            # try to award badges
            awards.instant(request)

            
            
            
            
            
            
            
            
            
            
