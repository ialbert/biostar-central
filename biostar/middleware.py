import datetime

MINIMUM_TIME = 60 * 1 # Time between updates to last visited

class LastVisitTimestampUpdaterMiddleware(object):
    def process_request(self, request):
        if request.user.is_authenticated():
            profile = request.user.get_profile()
            now = datetime.datetime.now()

            # Prevent writing to the database too often
            diff = (now - profile.last_visited).seconds
            if diff > MINIMUM_TIME:
                profile.last_visited = now
                profile.save()
    
        return None
