"""
Custom context processor
"""

from django.conf import settings

def extras(request):
    "Adds more data to each RequestContext"

    user = request.user
    
    # query- this is the primary query that goes into the search box
    q = request.REQUEST.get('q', '')
    
    # match - this is the secondary query value for matching other content
    m = request.REQUEST.get('m', '')
    

    return { 'BIOSTAR_VERSION': settings.BIOSTAR_VERSION, 'user':user, 'q':q, 'm':m }
