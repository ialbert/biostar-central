import datetime

class LastVisitTimpstampUpdaterMiddleware(object):
    def process_request(self, request):
        if request.user.is_authenticated():
            profile = request.user.get_profile()
            now = datetime.datetime.now()

            # Prevent writing to the database more often than 10 seconds
            diff = (now - profile.last_visited).seconds
            if diff > 10:
                profile.last_visited = now
                profile.save()
    
        return None
