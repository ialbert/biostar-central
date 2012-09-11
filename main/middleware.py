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
    SESSION_KEY, SORT_KEY, COUNT_KEY, TAB = "session-data", 'sortz', 'count', 'tabz'
    def __init__(self, request):
        self.request = request
        self.has_storage = request.user.is_authenticated()
        default = { self.COUNT_KEY:{ }, self.SORT_KEY:"rank", self.TAB:"all"}
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
        return self.data[self.TAB]
    
    def sort_order(self):
        "Stores the last sort order in the session"
        value = self.request.GET.get('sort', '').lower()
        value = value or self.data.get(self.SORT_KEY)
        if value not in SORT_MAP:
            value = "rank"
        self.data[self.SORT_KEY] = value
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
            
def generate_counts(request, weeks=6):
    "Returns the number of counts for each post type in the interval that has passed"
    user = request.user
    now  = datetime.now()
    
    key = 'countkey'
    if user.is_authenticated():
        since = user.profile.last_visited
    else:
        counts = cache.get(key)
        if counts:
            return counts
        since = now - timedelta(weeks=weeks)

    # posts since the last visit
    pairs = models.Post.objects.filter(type__in=POST_TOPLEVEL, status=POST_OPEN, creation_date__gt=since).order_by('-id').values_list("type", "answer_count")
    
    # establish how many of the posts have not been answered 
    values = [ p[0] for p in pairs ]
    unansw = len([ ptype for (ptype, pcount) in pairs if (ptype == POST_QUESTION) and (pcount == 0) ])
    
    # needs to be sorted for the groupby
    values.sort()

    # how many times does each post type appear in the list
    counts = dict( [ (POST_MAP[k], len(list(v))) for (k, v) in groupby(values) ] )
    
    # fill in unanswered posts
    counts['Unanswered'] = unansw
    
    if not user.is_authenticated():
        # store the cache key for non-authenticated users
        cache.set(key, counts, 600)

    return counts

class LastVisit(object):
    """
    Updates the last visit stamp at SESSION_UPDATE_TIME intervals
    """

    def process_request(self, request):
        user = request.user
        sess = Session(request)
        
        # check suspended status for users
        if user.is_authenticated() and (user.profile.status == USER_SUSPENDED):
            logout(request)
            messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
            return
            
        # handle anonymous users
        if not user.is_authenticated():
            # anonymous users
            request.user.can_moderate = False
            if request.path == "/":
                messages.info(request, 'Welcome to BioStar! Questions and Answers on Bioinformatics and Genomics!')
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
            
            # create nagging message for fixme posts
            fixme = models.Post.objects.filter(type=POST_FIXME, author=user, status=POST_OPEN)
            if fixme:
                first = fixme[0]
                messages.error(request, 'You have a post that does not conform the requirements. Please edit it: <a href="%s">%s</a>' % (first.get_absolute_url(), first.title)) 

            # remind user about voting every six weeks since the last vote
            # and only nag people with lower reputations ;-)
            if user.profile.score < 300:
                since = datetime.now() - timedelta(weeks=6)
                votes = models.Vote.objects.filter(author=user, date__gt=since)[:1]
                if not votes:
                    messages.info(request, '<i class="icon-info-sign"></i> Remember to <b>vote</b> on posts that you find useful!') 
                
            # try to award badges
            awards.instant(request)


class DebugMiddleware(object):
    """
    """

    def process_request(self, request):
        pass