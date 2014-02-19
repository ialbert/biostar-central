__author__ = 'ialbert'
from django.views.generic import DetailView, ListView, TemplateView, RedirectView, View
from haystack.views import SearchView
from haystack.forms import SearchForm
from haystack.query import SearchQuerySet
from haystack.utils import Highlighter

from django.conf import settings
from biostar.server.views import BaseListMixin
from ajax import ajax_error, ajax_success, ajax_error_wrapper, json_response

class SiteSearch(SearchView):
    extra_context = lambda x: dict(topic="search", page_title="Search")

class Search(BaseListMixin):
    template_name = "search/search.html"
    paginate_by = settings.PAGINATE_BY
    context_object_name = "results"
    page_title = "Search"

    def get_queryset(self):
        self.q = self.request.GET.get('q')
        query = SearchQuerySet().filter(content=self.q)

        return query

    def get_context_data(self, **kwargs):
        context = super(Search, self).get_context_data(**kwargs)
        context['q'] = self.q

        return context


def slow_highlight(query, text):
    "Invoked only if the search backend does not support highlighting"
    highlight = Highlighter(query)
    value = highlight.highlight(text)
    return value

#@ajax_error_wrapper
def search_title(request):
    "Handles title searches"
    q = request.GET.get('q','')

    results = SearchQuerySet().filter(content=q).highlight()[:25]

    items = []
    for row in results:
        ob = row.object
        snippet = '<br>'.join(x for x in row.highlighted['text'])
        snippet = snippet or slow_highlight(query=q, text=row.content)
        text = "%s" % row.title
        items.append(
            dict(id=ob.id, text=text, snippet=snippet, author=row.author),
        )

    msg = "OK"
    payload = dict(items=items)
    return json_response(payload)
    #return ajax_success(msg, items=items)