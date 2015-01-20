# Haystack search indices.
from biostar3.forum.models import Post
from django.db.models import Q
from haystack import indexes

# Create the search indices.
class PostIndex(indexes.SearchIndex, indexes.Indexable):
    text = indexes.CharField(document=True, use_template=True)
    title = indexes.CharField(model_attr='title')
    type = indexes.CharField(model_attr='type')
    content = indexes.CharField(model_attr='content')
    author = indexes.CharField(model_attr='author__name')
    vote_count = indexes.IntegerField(model_attr='vote_count')

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