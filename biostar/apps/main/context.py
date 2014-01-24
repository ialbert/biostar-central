__author__ = 'ialbert'
from django.conf import settings
from biostar import VERSION

def shortcuts(request):
    # These values will be added to each context
    context = {
        "GOOGLE_TRACKER": settings.GOOGLE_TRACKER,
        "BIOSTAR_VERSION": VERSION
    }
    return context