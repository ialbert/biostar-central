import logging
from django.conf import settings
from django_elasticsearch_dsl import fields, Index, Document
from django_elasticsearch_dsl_drf.compat import KeywordField, StringField
from django.conf import settings
from django_elasticsearch_dsl_drf.analyzers import edge_ngram_completion


class SpamDocument(Document):
    """Spam Elastic-search document."""

    class Index:
        name = settings.SPAM_INDEX_NAME

    id = fields.IntegerField(attr='id')
    is_spam = fields.BooleanField()

    title = StringField(
        fields={
            'raw': KeywordField(),
            'suggest': fields.CompletionField(),
            'edge_ngram_completion': StringField(
                analyzer=edge_ngram_completion
            ),
            'mlt': StringField(),
        }
    )

    content = StringField(
        fields={
            'raw': KeywordField(),
            'mlt': StringField(),
        }
    )


class TrainSpam(SpamDocument):

    class Index:
        name = settings.TRAIN_SPAM_INDEX
