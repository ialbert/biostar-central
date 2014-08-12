from django.db import transaction

from rest_framework import serializers, permissions

from biostar.server.ajax import perform_vote
from .models import Vote, Post


class VoteSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer class for `Vote` model.
    """
    class Meta:
        model = Vote
        fields = ('id', 'url', 'type', 'date', 'author', 'post')
        read_only_fields = ('date', 'author')

    def save(self, **kwargs):
        with transaction.atomic():
            msg = perform_vote(post=self.object.post,
                               user=self.object.author,
                               vote_type=self.object.type)
            if 'added' in msg.lower():
                self.object = Vote.objects.get(post=self.object.post,
                                               author=self.object.author,
                                               type=self.object.type)
            return self.object


class PostSerializer(serializers.HyperlinkedModelSerializer):
    """
    Serializer class for `Post` model.
    """
    class Meta:
        model = Post
        fields = ('id', 'url', 'title', 'author', 'votes', 'html')
                #TODO choose what fields to display


## Custom permissions #############################################################################
class IsOwnerOrReadOnly(permissions.BasePermission):
    """
    Custom permission to only allow author of a vote to delete it.
    """
    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        # Write permissions are only allowed to the author of the vote.
        return obj.author == request.user