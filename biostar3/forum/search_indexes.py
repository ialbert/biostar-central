# Haystack search indices.
from __future__ import absolute_import, division, print_function, unicode_literals
from biostar3.forum.models import Post, FederatedContent, BlogPost
from django.db.models import Q
from haystack import indexes
import json, time
from . import html

# Create the search indices.
class PostIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='title')
    type = indexes.CharField(model_attr='type')
    date = indexes.CharField()
    author = indexes.CharField()
    content = indexes.CharField()
    vote_count = indexes.IntegerField(model_attr='vote_count')
    view_count = indexes.IntegerField()
    reply_count = indexes.IntegerField(model_attr='reply_count')
    subtype = indexes.CharField(model_attr='subtype')
    domain = indexes.CharField()
    url = indexes.CharField()

    def get_model(self):
        return Post

    def prepare(self, obj):
        data = super(PostIndex, self).prepare(obj)
        data['boost'] = 1.0
        data['url'] = obj.get_absolute_url()
        data['domain'] = obj.root.usergroup.domain
        data['type'] = obj.get_type_display()
        data['view_count'] = obj.root.view_count
        data['date'] = obj.creation_date.strftime("%B %d, %Y")
        data['author'] = obj.author.name
        data['content'] = html.strip_tags(obj.content)
        return data

    def index_queryset(self, using=None):
        """
        Used when the entire index for model is updated.
        """
        cond = Q(type=Post.COMMENT) | Q(status=Post.DELETED)
        return self.get_model().objects.all().exclude(cond).select_related('root', "usergroup", "author")

    def get_updated_field(self):
        return "lastedit_date"


# Create the search indices for federated content.
class FederatedContentIndex(indexes.SearchIndex, indexes.Indexable):
    FIELDS = "title type vote_count domain url content".split()

    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField()
    type = indexes.CharField()
    url = indexes.CharField()
    date = indexes.CharField()
    content = indexes.CharField()
    vote_count = indexes.IntegerField()
    domain = indexes.CharField()
    author = indexes.CharField()

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

# Create the search indices for federated content.
class BlogPostIndex(indexes.SearchIndex, indexes.Indexable):
    FIELDS = "title type vote_count domain url content".split()

    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='title')
    url = indexes.CharField(model_attr='link')
    date = indexes.CharField(model_attr='creation_date')
    content = indexes.CharField()
    domain = indexes.CharField()
    type = indexes.CharField()
    author = indexes.CharField()

    def get_model(self):
        return BlogPost

    def prepare(self, obj):
        data = super(BlogPostIndex, self).prepare(obj)
        data['domain'] = obj.blog.usergroup.domain
        data['date'] = obj.creation_date.strftime("%B %d, %Y")
        data['content'] = html.strip_tags(obj.content)
        data['type'] = "Blog"
        data['author'] = obj.blog.title
        return data

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        query = self.get_model().objects.all().select_related("blog")
        return query

    def get_updated_field(self):
        return "date"