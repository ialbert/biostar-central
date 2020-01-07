from django_elasticsearch_dsl import Document
from django_elasticsearch_dsl.registries import registry
from .models import Post


class UserDocument(Document):
    pass


@registry.register_document
class PostDocument(Document):
    class Index:
        # Name of the Elasticsearch index
        name = 'post'
        # See Elasticsearch Indices API reference for available settings
        settings = {'number_of_shards': 1,
                    'number_of_replicas': 0}

    class Django:
        model = Post  # The model associated with this Document

        # The fields of the model you want to be indexed in Elasticsearch
        fields = [

            "status", "type", "title", "rank", "indexed", "is_toplevel", "answer_count", "accept_count",
            "reply_count", "comment_count", "vote_count", "thread_votecount", "view_count",
            "book_count", "subs_count", "creation_date", "lastedit_date", "sticky",
            "content", "html", "tag_val", "uid", "spam"
        ]

        # Ignore auto updating of Elasticsearch when a model is saved
        # or deleted:
        # ignore_signals = True

        # Don't perform an index refresh after every update (overrides global setting):
        # auto_refresh = False

        # Paginate the django queryset used to populate the index with the specified size
        # (by default it uses the database driver's default setting)
        # queryset_pagination = 5000
