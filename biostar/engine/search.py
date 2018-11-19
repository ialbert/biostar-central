import copy
import logging
from django.conf import settings
from whoosh.fields import Schema, ID, TEXT, STORED, NGRAM
from whoosh import index
from whoosh.qparser import MultifieldParser
from biostar.engine.models import Project, Job, Analysis, Data
from biostar.engine.const import *

logger = logging.getLogger('engine')


def get_obj_type(obj):

    type_string = ""

    if isinstance(obj, Job):
        type_string = "Job"
    elif isinstance(obj, Analysis):
        type_string = "Analysis"
    elif isinstance(obj, Data):
        type_string = "Data"
    elif isinstance(obj, Project):
        type_string = "Project"

    return type_string


def write_to_index(queryset=[]):
    """Create or update the index with objects in queryset. """

    schema = Schema(uid=ID(unique=True, stored=True), content=NGRAM(stored=True),
                    name=NGRAM(stored=True), type=STORED, url=STORED)

    # Index may already exists in the root directory
    if index.exists_in(settings.INDEX_ROOT):
        ix = index.open_dir(dirname=settings.INDEX_ROOT)
    else:
        ix = index.create_in(dirname=settings.INDEX_ROOT, schema=schema)

    writer = ix.writer()
    for obj in queryset:
        # update_document updates if it exists otherwise creates new 'document'
        if obj.project.is_public:
            writer.update_document(uid=obj.uid,
                                   name=obj.name,
                                   content=obj.text,
                                   type=get_obj_type(obj=obj),
                                   url=obj.url())
    writer.commit()

    return writer


def search_index(text_query):
    """Search index for a given text."""

    # Index has not bee
    if not index.exists_in(settings.INDEX_ROOT):
        logger.error("Index has not been built")
        return

    ix = index.open_dir(dirname=settings.INDEX_ROOT)

    # Parse the user query, looking to match the "name" or "content"
    # fields of all indexed objects.
    parser = MultifieldParser(["name", "content"], schema=ix.schema)
    query = parser.parse(text_query)

    results = []
    with ix.searcher() as searcher:
        all_results = searcher.search(query, limit=None)
        all_results.fragmenter.surround = MAX_NAME_LEN

        for res in all_results:
            res_copy = copy.deepcopy(dict(res.items()))
            res_copy["name"] = res.highlights("name") or res_copy["name"]
            res_copy["content"] = res.highlights("content") or res_copy["content"]
            results.append(res_copy)

    return results









