import datetime
from django.contrib.auth import logout
from django.contrib import messages
from main.server import models, notegen, html
from main.server.const import *

class LastVisit(object):
    """
    Updates the last visit stamp at MINIMUM_TIME intervals
    """
    # minimum elapsed time
    MINIMUM_TIME = 60 * 10 # every 5 minutes

    def process_request(self, request):
        
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
            if diff > self.MINIMUM_TIME:
               
                last = user.profile.last_visited
               
                questions = models.Post.objects.filter(type=POST_QUESTION, creation_date__gt=last).count()
                unanswered = models.Post.objects.filter(type=POST_QUESTION, answer_count=0,  creation_date__gt=last).count()
                guides = models.Post.objects.filter(type=POST_GUIDE, answer_count=0, creation_date__gt=last).count()
                planet = models.Post.objects.filter(type=POST_BLOG,  creation_date__gt=last).count()
                
                counts = dict(planet=planet,  unanswered=unanswered, questions=questions, guides=guides)
                request.session[SESSION_POST_COUNT] = counts
                
                models.UserProfile.objects.filter(user=user).update(last_visited=now)
               
                # award the beta tester badge
                #models.apply_award(request=request, user=user, badge_name=BETA_TESTER_BADGE, messages=messages)
                
                # recompute the posts
                
            # a handy shortcut
            request.user.can_moderate = profile.can_moderate
        else:
            request.user.can_moderate = False

        return None


class PermissionsMiddleware(object):
    ''' Calculates the logged-in user's permissions and adds it to the request object. '''
    def process_request(self, request):
        pass