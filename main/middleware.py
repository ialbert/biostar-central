import datetime

class LastVisit(object):
    """
    Updates the last visit stamp at MINIMUM_TIME intervals
    """
    # minimum elapsed time
    MINIMUM_TIME = 60 * 1 # 1 minute

    def process_request(self, request):
        if request.user.is_authenticated():
            profile = request.user.get_profile()
            now = datetime.datetime.now()
            diff = (now - profile.last_visited).seconds
            
            # Prevent writing to the database too often
            if diff > self.MINIMUM_TIME:
                profile.last_visited = now
                profile.save()
    
        return None


class PermissionsMiddleware(object):
    ''' Calculates the logged-in user's permissions and adds it to the request object. '''
    def process_request(self, request):
        if request.user.is_authenticated():
            request.permissions = request.user.profile.permissions
        else:
            request.permissions = []