from __future__ import absolute_import, division, print_function, unicode_literals

from django.conf import settings
from django.contrib import messages
from django.core.cache import cache
from biostar3 import VERSION

def shortcuts(request):
    # These values will be added to each context

    context = {
        "BIOSTAR_VERSION": VERSION,
    }

    return context

def modify_context(ctx, request):
    """
    Mutates the context and populates sort, limit and query fields.
    """
    sort = request.GET.get('sort', '')
    limit = request.GET.get('limit', '')
    q = request.GET.get('q', '')

    if sort and sort not in settings.POST_SORT_MAP:
        messages.warning(request, settings.POST_SORT_INVALID_MSG)
        sort = ''

    ctx['sort'] = sort
    ctx['limit'] = limit
    ctx['q'] = q[:100]
