__author__ = 'ialbert'
from biostar3.utils.compat import *
from haystack.query import SearchQuerySet, AutoQuery
from haystack.utils import Highlighter
import logging

logger = logging.getLogger(__name__)

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


def plain(query, content_type=None):
    """
    Simplest search query.
    """
    content = AutoQuery(query)
    results = SearchQuerySet().filter(content=content).highlight()[:100]

    if content_type:
        results = filter(lambda x: x.content_type() == content_type, results)

    for row in results:
        try:
            if not row:
                continue
            context = join_highlights(row)
            context = context or slow_highlight(query=query, text=row.content)
            row.context = context
        except Exception as exc:
            # Occasionally haystack triggers errors. row.obj is None - perhaps for deleted objects
            # that are still in the search index.
            logger.error(exc)

    return results




