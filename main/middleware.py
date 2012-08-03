import datetime
from itertools import *
from django.contrib.auth import logout
from django.contrib import messages
from django.conf import settings
from main.server import models, notegen, html, awards
from main.server.const import *

settings.CONTENT_INDEXING = True

class Session(object):
    """
    This object maintains various session for user. Also works for non-authenticated users and avoids creating database sessions for them
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
        
    def counts(self):
        return self.data[self.COUNT_KEY]
     
    def tabpill(self, value=None):
        "Facilitates navigation by remebering the last visited tab and pill"
        
        # these are the old values
        otab, opill = self.data[self.TAB], self.data[self.PILL]
     
        # nothing coming in, keep the old values
        if not value:
            return(otab, opill)
        
        # the tab is always set, the pill is only set
        # on valid requests
        self.data[self.TAB] = value
        
        # hitting the post, select last pill
        if value == "posts":
            tab, pill = value, opill
            
        # a valid tab other than posts
        elif value in VALID_TABS:
            return (value, "")
            
        # request for a valid pill link
        elif value in VALID_PILLS:
            tab, pill = "posts", value
        
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

def get_counts(since):
    "Returns the number of counts for each post type in the interval that has passed"
    
    # get the posts with a sanity limit
    values = models.Post.objects.filter(type__in=POST_TOPLEVEL, creation_date__gt=since).values_list("type", flat=True)[:1000]
    
    # how many times does each post type appear in the list
    counts = dict( [ (POST_MAP[k], len(list(v))) for (k, v) in groupby(values) ] )
    
    return counts

class LastVisit(object):
    """
    Updates the last visit stamp at MINIMUM_TIME intervals
    """

    def process_request(self, request):
        
        # content generated within the last three months
        
        now = datetime.datetime.now() 
        since = now - datetime.timedelta(weeks=50)
        
        counts = get_counts(since)
                
        if request.user.is_authenticated():
            user = request.user
            profile = user.get_profile()
            
            if profile.status == USER_SUSPENDED:
                logout(request)
                messages.error(request, 'Sorry, this account has been suspended. Please contact the administrators.')
                return None
            
            now = datetime.datetime.now()
            diff = (now - profile.last_visited).seconds
            
            # Prevent writing to the database too often
            if diff > settings.SESSION_UPDATE_TIME:
               
               # create nagging message for fixme posts
                fixme = models.Post.objects.filter(type=POST_FIXME, author=user, status=POST_OPEN)
                if fixme:
                    first = fixme[0]
                    messages.error(request, 'You have a post that does not conform the requirements. Please edit it: <a href="%s">%s</a>' % (first.get_absolute_url(), first.title) )
            
                last = user.profile.last_visited
               
                pairs = [
                    ('questions', POST_QUESTION), ('forum', POST_FORUM), ('tutorials', POST_TUTORIAL),
                    ('planet', POST_BLOG), ('jobs', POST_JOB), ('tools', POST_TOOL)
                ]
                
                counts = {}
                for key, value in pairs:
                    counts[key] = models.Post.objects.filter(type=value, creation_date__gt=last).count()
                
                # this is a special case
                counts['unanswered'] = models.Post.objects.filter(type=POST_QUESTION, status=POST_OPEN, answer_count=0,  creation_date__gt=last).count()
                
                request.session[SESSION_POST_COUNT] = counts 
                
                user.profile.last_visited = now
                user.profile.save()
               
                awards.instant(request)
                
                # add the visit to the database
                try:
                    # trying to establish the IP location
                    ip1 = request.META.get('REMOTE_ADDR', '')
                    ip2 = request.META.get('HTTP_X_FORWARDED_FOR','').split(",")[0].strip()
                    ip  = ip1 or ip2 or '0.0.0.0'
                    models.Visit.objects.create(ip=ip, user=user)
                except Exception ,exc:
                    print '*** ip handling error %s' % exc
                    
            # a handy shortcut
            request.user.can_moderate = profile.can_moderate
            
        else:
            request.user.can_moderate = False
            has_session = request.session.get(FIRST_SESSION)
            if not has_session:
                request.session[FIRST_SESSION] = True
                messages.info(request, 'Welcome to BioStar! Questions and answers on bioinformatics, computational genomics and systems biology.')

        return None

class ErrorCheckMiddleware(object):
    ''' Calculates the logged-in user's permissions and adds it to the request object. '''
    
    def process_exception(self, request, exc):
        path = request.path
        params = html.Params(exc=exc, path=path)
        return html.template(request, name='500.html', params=params)
