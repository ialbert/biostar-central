__author__ = 'ialbert'
from biostar.apps.posts.models import Post
from biostar.apps.planet.models import BlogPost

from haystack import indexes

# Create the search indices.
class PostIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='title')
    type = indexes.CharField(model_attr='type')
    content = indexes.CharField(model_attr='content')
    author = indexes.CharField(model_attr='author__name')

    def get_model(self):
        return Post

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return self.get_model().objects.exclude(status=Post.DELETED)

    def get_updated_field(self):
        return "lastedit_date"

# Create the search indices.
class BlogPostIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='get_title')
    content = indexes.CharField(model_attr='html')
    author = indexes.CharField(model_attr='blog__title')

    def get_model(self):
        return BlogPost

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        query = self.get_model().objects.all()
        return query

    def get_updated_field(self):
        return "creation_date"
