"""
Custom context processor
"""
import itertools

from django.conf import settings

from main.server import models

def extras(request):
    "Adds more data to each RequestContext"

    user = request.user
    
    # query- this is the primary query that goes into the search box
    q = request.REQUEST.get('q', '')
    
    # match - this is the secondary query value for matching other content
    m = request.REQUEST.get('m', '')
    

    return { 'BIOSTAR_VERSION': settings.BIOSTAR_VERSION, 'user':user, 'q':q, 'm':m }

def popular_tags(request):
    tags1 = models.Tag.objects.filter(name='galaxy') # Special treatment for Galaxy folks
    tags2 = models.Tag.objects.exclude(name='galaxy').order_by('-count')[:5]
    tags = itertools.chain(tags1, tags2)

    return {'popular_tags' : tags}


