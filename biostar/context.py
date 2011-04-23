"""
Custom context processor
"""

from django.conf import settings 

def extras(request):
    "Adds more data to each RequestContext"
    return {'BIOSTAR_VERSION': settings.BIOSTAR_VERSION, 'permissions':request.permissions}
