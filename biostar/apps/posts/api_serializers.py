from django.db import transaction

from rest_framework import serializers, permissions

from biostar.server.ajax import perform_vote, validate_vote, VALID_VOTE, DOWNVOTE, \
UPVOTED_OWN_POST, ACCEPTED_NOT_ANSWER, ACCEPTED_NOT_OWN_QUESTION, VOTE_VALIDATION_MSGS
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
        """
        Override the original `save` method in order to use `perform_vote` when storing a new
        vote. `perform_vote` takes care of the side effects.
        """
        with transaction.atomic():
            msg = perform_vote(post=self.object.post,
                               user=self.object.author,
                               vote_type=self.object.type)
            if 'added' in msg.lower():
                self.object = Vote.objects.get(post=self.object.post,
                                               author=self.object.author,
                                               type=self.object.type)
            return self.object

    # Cached result of the validation.
    _validation_result = None

    def get_validation_result(self, attrs):
        """
        The actual validation of the vote is executed by `validate_vote`. This method calls once
        the actual validation and caches the result in `_validation_result`.
        """
        if self._validation_result == None:
            request = self.context['request']
            vote = Vote(author=request.user, post=attrs['post'], type=int(attrs['type']))
            self._validation_result = validate_vote(vote)
        return self._validation_result

    def validate(self, attrs):
        """
        Object-level validation of the `Vote` object. This method is executed if all the
        field-level validation methods (like `validate_type()`) have passed.
        """
        # A reminder: from a serializer you can access the relative vie and its request using
        # the context attribute.
        # self.context = {
        #     u'view': <biostar.apps.posts.views.VoteViewSet object at 0x10e554090>,
        #     u'request': <rest_framework.request.Request object at 0x10f2848d0>,
        #     u'format': None
        # }
        # request = self.context['request']

        validation_result = self.get_validation_result(attrs)
        if not validation_result == VALID_VOTE:
            raise serializers.ValidationError(VOTE_VALIDATION_MSGS[validation_result] + '.')
        return attrs

    def validate_type(self, attrs, source):
        """
        Validation of `type` attribute.
        """
        validation_result = self.get_validation_result(attrs)

        # The only validation error related to this field is downvote.
        if validation_result == DOWNVOTE:
            raise serializers.ValidationError(VOTE_VALIDATION_MSGS[validation_result] + '.')
        return attrs

    def validate_post(self, attrs, source):
        """
        Validation of `post` attribute.
        """
        validation_result = self.get_validation_result(attrs)

        # The validation errors related to this field are: UPVOTED_OWN_POST, ACCEPTED_NOT_ANSWER,
        # ACCEPTED_NOT_OWN_QUESTION.
        if validation_result in (UPVOTED_OWN_POST, ACCEPTED_NOT_ANSWER, ACCEPTED_NOT_OWN_QUESTION):
            raise serializers.ValidationError(VOTE_VALIDATION_MSGS[validation_result] + '.')
        return attrs


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