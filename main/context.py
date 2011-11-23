"""
Custom context processor
"""

from django.conf import settings

def extras(request):
    "Adds more data to each RequestContext"

    user = request.user
         
    return { 'BIOSTAR_VERSION': settings.BIOSTAR_VERSION, 'user':user }
