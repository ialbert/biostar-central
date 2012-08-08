from datetime import datetime, timedelta
from itertools import *
from django.contrib.auth import logout
from django.contrib import messages
from django.conf import settings
from main.server import models, notegen, html, awards
from main.server.const import *

settings.CONTENT_INDEXING = True

class Session(object):
    """
    This object maintains various session data for a user. To minimize database access
    saves all at just one time. It also works for non-authenticated users but avoids
    creating a database sessions for them
    """
    SESSION_KEY, SORT_KEY, COUNT_KEY, TAB, PILL = "session-data", 'sort', 'count', 'tab', 'pill'
    def __init__(self, request):
        self.request = request
        self.has_storage = request.user.is_authenticated()
        default = { self.COUNT_KEY:{ }, self.SORT_KEY:"rank", self.TAB:"posts", self.PILL:"all" }
        if self.has_storage:
            self.data = self.request.session.get(self.SESSION_KEY, default )
        else:
            self.data = default
    
    def save(self):
        "Saves the counts back to the session"
        if self.has_storage:
            self.request.session[self.SESSION_KEY] = self.data
             
    def tabpill(self, value=None):
        "Facilitates navigation by remebering the last visited tab and pill"
        
        # these are the old values
        otab, opill = self.data[self.TAB], self.data[self.PILL]
     
        # nothing specified, keep the old values
        if not value:
            return(otab, opill)
        
        # the tab must always be set
        self.data[self.TAB] = value
        
        # requesting the post tab, select and return the last pill
        if value == "posts":
            return (value, opill)
            
        # a valid tab other than posts
        elif value in VALID_TABS:
            return (value, "")
            
        # request for a valid pill link
        elif value in VALID_PILLS:
            tab, pill = "posts", value
        
        # only set the pill when a valid PILL request comes in
        self.data[self.TAB]  = tab
        self.data[self.PILL] = pill
       
        return tab, pill
    
    def sort_order(self):
        "Stores the last sort order in the session"
        value = self.request.GET.get('sort', '').lower()
        last  = self.data.get(self.SORT_KEY, '')
        value = value or last
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
            
def generate_counts(request, weeks=50):
    "Returns the number of counts for each post type in the interval that has passed"
    user = request.user
    now  = datetime.now()
    
    if user.is_authenticated():
        since = user.profile.last_visited
    else:
        since = now - timedelta(weeks=weeks)
    
    # the the posts since the last time
    values = models.Post.objects.filter(type__in=POST_TOPLEVEL, creation_date__gt=since).values_list("type", flat=True)[:1000]
        
    # how many times does each post type appear in the list
    counts = dict( [ (POST_MAP[k], len(list(v))) for (k, v) in groupby(values) ] )
    
    # fill in unanswered posts
    counts['Unanswered'] = models.Post.objects.filter(type=POST_QUESTION, status=POST_OPEN, answer_count=0,  creation_date__gt=since).count()
     
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

            # try to award badges
            awards.instant(request)

