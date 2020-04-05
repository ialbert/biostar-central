import logging
from django.conf import settings
from django_elasticsearch_dsl import fields, Index, Document
from django.conf import settings


class SpamDocument(Document):
    """Spam Elastic-search document."""

    def __init__(self, indexname=None, *args, **kwargs):
        indexname = indexname or settings.SPAM_INDEX_NAME
        super(SpamDocument, self).__init__(*args, **kwargs)
        self.Index.name = indexname

    class Index:
        name = settings.SPAM_INDEX_NAME

    uid = fields.TextField(attr='id')
    title = fields.TextField()
    content = fields.TextField()
    is_spam = fields.BooleanField()
