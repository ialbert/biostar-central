
import logging
from functools import reduce
from django.db.models import Q
from django.utils.text import smart_split
from django.conf import settings
from django import forms

from biostar.recipes.auth import get_project_list
from biostar.recipes.models import Job, Analysis, Data
from biostar.recipes.const import *

logger = logging.getLogger('engine')


def search(request):
    """
    Searches recipes
    """

    results = []

    # Search each model type by title, text, owner email/name.
    search_fields = ['name', 'text', 'owner__email', 'owner__profile__name']

    # Get the objects the user can access
    projects = get_project_list(user=request.user)

    recipes = Analysis.objects.filter(project__in=projects, root=None, deleted=False)
    # Load query from GET request.
    search_form = SearchForm(queryset=recipes, search_fields=search_fields,
                             data=request.GET or {})

    # Add search results to dict
    if search_form.is_valid():
        results = search_form.get_queryset()

    return results


class SearchForm(forms.Form):

    def __init__(self, queryset=None, search_fields=None, *args, **kwargs):
        self.queryset = queryset
        self.search_fields = search_fields
        super(SearchForm, self).__init__(*args, **kwargs)

    q = forms.CharField(label='Search', required=False)

    def clean_q(self):
        query = self.cleaned_data['q'].strip()

        if len(query) <= settings.SEARCH_CHAR_MIN:
            raise forms.ValidationError("Enter more than 3 characters.")

        return query

    def get_queryset(self):
        qs = self.queryset
        query = self.cleaned_data.get('q')

        if query:
            qs = qs.filter(search_filter(self.search_fields, query))

        return qs


def search_filter(search_fields, query_string):
    """search_fields example: ['name', 'category__name', '@description', '=id']
    """

    query_string = query_string.strip()

    filters = []
    first = True

    for bit in split_text_query(query_string):

        queries = [Q(**{search_param(field_name, first): bit}) for field_name in search_fields]
        filters.append(reduce(Q.__or__, queries))
        first = False

    return reduce(Q.__and__, filters) if len(filters) else Q(pk=None)


def search_param(field_name, is_first_word):
    if field_name.startswith('^') and is_first_word:
        return "%s__istartswith" % field_name[1:]
    elif field_name.startswith('@'):
        return "%s__search" % field_name[1:]
    elif field_name.startswith('='):
        return "%s__iexact" % field_name[1:]
    else:
        return "%s__icontains" % field_name


def split_text_query(query):
    """Filter out stopwords but only if there are useful words"""

    split_query = list(smart_split(query))
    filtered_query = [bit for bit in split_query if bit not in STOPWORDS]

    return filtered_query if len(filtered_query) else split_query



