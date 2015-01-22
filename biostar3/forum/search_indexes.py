# Haystack search indices.
from biostar3.forum.models import Post, FederatedContent
from django.db.models import Q
from haystack import indexes
import json

# Create the search indices.
class PostIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='title')
    type = indexes.CharField(model_attr='type')
    content = indexes.CharField(model_attr='content')
    vote_count = indexes.IntegerField(model_attr='vote_count')
    domain = indexes.CharField()
    url = indexes.CharField()

    def prepare_url(self, obj):
        return "url"

    def prepare_domain(self, obj):
        return "local"

    def prepare_type(self, obj):
        return "%s" % obj.get_type_display()

    def get_model(self):
        return Post

    def prepare(self, obj):
        data = super(PostIndex, self).prepare(obj)
        data['boost'] = 1.0
        return data

    def index_queryset(self, using=None):
        """
        Used when the entire index for model is updated.
        """
        cond = Q(type=Post.COMMENT) | Q(status=Post.DELETED)
        return self.get_model().objects.all().exclude(cond)

    def get_updated_field(self):
        return "lastedit_date"

# Create the search indices for federated content.
class FederatedContentIndex(indexes.SearchIndex, indexes.Indexable):
    FIELDS = "title type vote_count domain url content".split()

    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField()
    type = indexes.CharField()
    url = indexes.CharField()
    content = indexes.CharField()
    vote_count = indexes.IntegerField()
    domain = indexes.CharField()

    def prepare(self, obj):
        self.prepared_data = super(FederatedContentIndex, self).prepare(obj)

        fields = json.loads(obj.content)

        for field in self.FIELDS:
            self.prepared_data[field] = fields[field]

        return self.prepared_data

    def get_model(self):
        return FederatedContent

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        query = self.get_model().objects.all()
        return query

    def get_updated_field(self):
        return "changed"