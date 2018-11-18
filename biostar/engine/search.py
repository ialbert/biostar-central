
from django.conf import settings
from whoosh.fields import Schema, ID, TEXT, STORED
from whoosh import index
from whoosh.qparser import MultifieldParser
from biostar.engine.models import Project, Job, Analysis, Data


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


def get_schema():
    """Return universally used search index schema """

    return Schema(uid=ID(unique=True, stored=True), content=TEXT(stored=True),
                  name=TEXT(stored=True), type=STORED)


def write_to_index(queryset=[]):
    """Create or update the index with objects in queryset. """

    # Index may already exists in the root directory
    if index.exists_in(settings.INDEX_ROOT):
        ix = index.open_dir(dirname=settings.INDEX_ROOT)
    else:
        ix = index.create_in(dirname=settings.INDEX_ROOT, schema=get_schema())

    writer = ix.writer()
    for obj in queryset:
        # update_document updates if it exists otherwise creates new 'document'
        writer.update_document(uid=obj.uid,
                               name=obj.name,
                               content=obj.text,
                               type=get_obj_type(obj=obj))
    writer.commit()

    return writer


def search_index(text_query):

    ix = index.open_dir(dirname=settings.INDEX_ROOT)

    # Parse the user query
    parser = MultifieldParser(["name", "content"], schema=ix.schema)

    print(text_query, "TEXT")
    query = parser.parse(u"{}".format(text_query))

    with ix.searcher() as searcher:
        results = searcher.search(query, limit=None)

    print([x for x in ix.searcher().documents()])
    print(results, ix, len(results))

    return results









