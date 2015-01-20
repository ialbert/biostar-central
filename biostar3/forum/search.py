__author__ = 'ialbert'
from haystack.query import SearchQuerySet, AutoQuery
from haystack.utils import Highlighter


def slow_highlight(query, text):
    "Invoked only if the search backend does not support highlighting"
    highlight = Highlighter(query)
    value = highlight.highlight(text)
    return value


def join_highlights(row):
    "Joins the highlighted text"
    if type(row.highlighted) is dict:
        return ''

    # Unable to highlight by the back end
    if not row.highlighted:
        return ''

    return '<br>'.join(x for x in row.highlighted)


def plain(query):
    content = AutoQuery(query)
    results = SearchQuerySet().filter(content=content).highlight()[:100]

    for row in query:
        #context = join_highlights(row)
        #context = context or slow_highlight(query=query, text=row.content)
        #row.context = ''
        pass

    return results




