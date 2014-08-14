import copy

from django.contrib.sites.models import get_current_site
from django.db import transaction
from django.core import exceptions

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

    def save_object(self, obj, **kwargs):
        """
        Override the original `save_object` method in order to use `perform_vote` when saving a new
        vote. `perform_vote` takes care of the side effects like updating the user's reputation.

        `obj` refers to `self.object` which is the object that was just saved (and displayed in
        the response).
        """
        with transaction.atomic():
            msg = perform_vote(post=obj.post,
                               user=obj.author,
                               vote_type=obj.type)
            # We should have done something like:
            # obj.save()
            # So the object is saved to the db and `self.object` points to it.
            # But since we have to use `perform_vote` to save a vote, then we have to manually
            # fetch the vote from the db and assign it to `self.object`.
            #
            # `perform_vote`, among the other things, implements also the following rule:
            # Rule V5: if a vote with the same author, post, type does exist in the database, the
            # vote itself is deleted.

            if 'added' in msg.lower():
                self.object = Vote.objects.get(post=obj.post,
                                               author=obj.author,
                                               type=obj.type)

    # Cached result of the validation.
    _validation_result = None

    def get_validation_result(self, attrs):
        """
        The actual validation of the vote is executed by `validate_vote`. This method calls once
        the actual validation and caches the result in `_validation_result`.
        """
        # Rule V1: downvotes not allowed.
        # Rule V2: a user can not upvote her own post.
        # Rule V3: accept votes are only allowed for posts of type "answer".
        # Rule V4: the author of an accept vote must match the author of the root post (the
        # original question).
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
    permalink = serializers.SerializerMethodField('get_absolute_url')
    top_level = serializers.BooleanField(source='is_toplevel', read_only=True)

    class Meta:
        model = Post
        fields = ('id', 'url', 'title', 'author', 'lastedit_user', 'rank', 'status', 'type',
                  'top_level', 'vote_count', 'votes', 'view_count', 'reply_count', 'comment_count',
                  'book_count', 'subs_count', 'thread_score', 'creation_date', 'lastedit_date',
                  'sticky', 'has_accepted', 'root', 'parent', 'children', 'tag_val', 'content',
                  'html', 'permalink',)
        read_only_fields = ('author', 'lastedit_user', 'rank', 'status', 'vote_count', 'votes',
                            'view_count', 'reply_count', 'comment_count', 'book_count',
                            'subs_count', 'thread_score', 'creation_date', 'lastedit_date',
                            'sticky', 'has_accepted', 'root', 'children', 'html',)

    def get_absolute_url(self, obj):
        request = self.context['request']
        return 'http://{}{}'.format(get_current_site(request).domain, obj.get_absolute_url()),

    def validate(self, attrs):
        """
        Object-level validation of the `Post` object. This method is executed if all the
        field-level validation methods (like `validate_title()`) have passed.
        """
        post = self.get_new_post_preview(attrs)

        # Rule P6: answers allowed only if the parent post is open.
        if (post.type == Post.ANSWER and
            post.parent.status != Post.OPEN):
            raise serializers.ValidationError("Only open posts can have answers.")

        # Rule P7: comments allowed only if the root post is open or closed.
        if (post.type == Post.COMMENT and
            post.parent.root.status not in [Post.OPEN, Post.CLOSED]):
            raise serializers.ValidationError("Only open and closed posts can have comments.")

        return attrs

    def validate_title(self, attrs, source):
        """
        Validation of `title` attribute.
        """
        # Rule P1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        # Rule P2: let x be the title length: 10 <= x <= 200; x >= 3 words.
        from .views import valid_title
        try:
            valid_title(attrs[source])
        except exceptions.ValidationError as ex:
            raise serializers.ValidationError(ex.message)
        if len(attrs[source]) > 200:
            raise serializers.ValidationError("The tile must be shorter than 200 chars.")

        return attrs

    def validate_type(self, attrs, source):
        """
        Validation of `type` attribute.
        """
        # Rule P1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        return attrs

    def validate_parent(self, attrs, source):
        """
        Validation of `parent` attribute.
        """
        # Rule P1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        # Rule P3: consistent parents:
        #  - The parent of a toplevel post must be emtpy or the post itself;
        #  - There must be a parent for non toplevel posts;
        #  - The parent of an answer must be a toplevel post.
        self.rule_consistent_parent(attrs, source)

        return attrs

    def validate_tag_val(self, attrs, source):
        """
        Validation of `tag_val` attribute.
        """
        # Rule P4: let x be the tag_val length: 0 < x <= 50, x <= 5 words.
        from .views import valid_tag
        try:
            valid_tag(attrs[source])
        except exceptions.ValidationError as ex:
            raise serializers.ValidationError(ex.message)

        return attrs

    def validate_content(self, attrs, source):
        """
        Validation of `content` attribute.
        """
        # Rule P5: let x be the content length: 80 < x < 15000.
        if len(attrs[source]) < 80 or len(attrs[source]) > 15000:
            raise serializers.ValidationError("The content must be longer than 80 chars "
                                              "and shorter than 15000 chars.")

        return attrs

    def get_new_post_preview(self, attrs):
        post = copy.copy(self.object) or Post()

        for key, val in attrs.items():
            setattr(post, key, val)

        # Fix the type: the filed `type` in `attrs` is a string like '0', but it must be an
        # integer in the real Post instance.
        post.type = int(post.type)

        return post

    def rule_field_not_updatable(self, attrs, source):
        post = self.get_new_post_preview(attrs)

        if self.context['request'].method.upper() in ['PUT', 'PATCH']:
            if getattr(post, source) != getattr(self.object, source):
                raise serializers.ValidationError("You can not update this field.")

    def rule_consistent_parent(self, attrs, source):
        post = self.get_new_post_preview(attrs)

        # The parent of a toplevel post must be emtpy or the post itself.
        if post.is_toplevel:
            if post.parent and not post.parent == post:
                raise serializers.ValidationError("The parent of a top level post must be either"
                                                  " emtpy or the post itself.")
        # There must be a parent for non toplevel posts.
        else:
            if not post.parent:
                raise serializers.ValidationError("Non top level posts must have parents.")

        # The parent of an answer must be a toplevel post.
        if post.type == Post.ANSWER:
            if not post.parent.is_toplevel:
                raise serializers.ValidationError("Answers must have top level parents.")


## Custom permissions #############################################################################
class IsOwnerOrReadOnly(permissions.BasePermission):
    """
    Custom permission to allow only the author of a vote to delete it.
    """
    # Rule V6: only the author of a vote can delete it.
    # Rule P8: users can update/delete only posts they authored.
    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        # Write permissions are only allowed to the author of the vote.
        return obj.author == request.user


class IsOpenOrReadOnly(permissions.BasePermission):
    """
    Custom permission to allow updates only for open posts.
    """
    # Rule P9: only open posts can be updated.
    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        # Write permissions are only allowed to the open posts.
        return obj.status == Post.OPEN