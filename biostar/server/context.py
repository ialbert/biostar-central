__author__ = 'ialbert'
from django.conf import settings
from biostar import VERSION

def shortcuts(request):
    # These values will be added to each context
    context = {
        "GOOGLE_TRACKER": settings.GOOGLE_TRACKER,
        "SITE_STYLE_CSS": settings.SITE_STYLE_CSS,
        "TOPICS": settings.DEFAULT_TOPICS,
        "BIOSTAR_VERSION": VERSION,
        "TRAFFIC": 0,
    }
    return context