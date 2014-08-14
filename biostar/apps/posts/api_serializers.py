import copy

from django.contrib.sites.models import get_current_site
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

        # Rule x: answers allowed only if the parent post is open.
        if (post.type == Post.ANSWER and
            post.parent.status != Post.OPEN):
            raise serializers.ValidationError("Only open posts can have answers.")

        # Rule x: comments allowed only if the root post is open or closed.
        if (post.type == Post.COMMENT and
            post.parent.root.status not in [Post.OPEN, Post.CLOSED]):
            raise serializers.ValidationError("Only open and closed posts can have somments.")

        return attrs

    def validate_title(self, attrs, source):
        """
        Validation of `title` attribute.
        """
        # Rule 1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        #from .views import valid_title
        #
        #try:
        #    valid_title(attrs[source])
        #except exceptions.ValidationError as ex:
        #    raise serializers.ValidationError(ex.message)
        #
        ### If it is an update/partial_update request and the post is not is_toplevel then the title
        ### cannot be changed.
        ##if (self.context['request'].method.upper() in ['PUT', 'PATCH'] and
        ##    not self.object.is_toplevel and
        ##    self.object.title != attrs[source]):
        ##    raise serializers.ValidationError('You can edit titles only for top level posts.')

        return attrs

    def validate_type(self, attrs, source):
        """
        Validation of `type` attribute.
        """
        # The `type` field is numeric, but here it is serialized as a string (like '0').
        # I am creating a new `attrs_fixed` with a nuumeric `type` field.
        #attrs_fixed = attrs.copy()
        #attrs_fixed[source] = int(attrs_fixed[source])

        # Rule 1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        ## If it is an update/partial_update request and the post is not is_toplevel then the type
        ## cannot be changed.
        #if (self.context['request'].method.upper() in ['PUT', 'PATCH'] and
        #    not self.object.is_toplevel and
        #    self.object.type != int(attrs[source])):
        #    raise serializers.ValidationError('You can edit types only for top level posts.')

        return attrs

    def validate_parent(self, attrs, source):
        """
        Validation of `parent` attribute.
        """
        # Rule 1: users can only update the `content` and `tag_val` fields.
        self.rule_field_not_updatable(attrs, source)

        # Rule 2: consistent parents.
        self.rule_consistent_parent(attrs, source)

        ## If the post is a toplevel post, then its parent must be the post itself.
        #if (int(attrs['type']) in Post.TOP_LEVEL and
        #    not (attrs['parent'] is self.object or not attrs['parent'])):
        #    raise serializers.ValidationError('The parent of a top level post must be emtpy or '
        #                                      'the post itself.')
        ## If the post is not a toplevel post, there must be a parent and the parent must be a top
        ## level post if the current post is an answer.
        #if int(attrs['type']) not in Post.TOP_LEVEL:
        #    if not attrs['parent']:
        #        raise serializers.ValidationError('Non top level posts must have a parent.')
        #    if (int(attrs['type']) == Post.ANSWER and
        #        not attrs['parent'].type in Post.TOP_LEVEL):
        #        raise serializers.ValidationError('The parent of an answer must be a top level '
        #                                          'post.')
        #
        ## If it is an update/partial_update request and the post is not is_toplevel then the parent
        ## cannot be changed.
        #if (self.context['request'].method.upper() in ['PUT', 'PATCH'] and
        #    not self.object.is_toplevel and
        #    self.object.parent != attrs[source]):
        #    raise serializers.ValidationError('You can edit parents only for top level posts.')

        return attrs

    def validate_tag_val(self, attrs, source):
        """
        Validation of `tag_val` attribute.
        """
        #from .views import valid_tag
        #
        #try:
        #    valid_tag(attrs[source])
        #except exceptions.ValidationError as ex:
        #    raise serializers.ValidationError(ex.message)
        #
        ## If it is an update/partial_update request and the post is not is_toplevel then the
        ## tag_val cannot be changed.
        #if (self.context['request'].method.upper() in ['PUT', 'PATCH'] and
        #    not self.object.is_toplevel and
        #    self.object.tag_val != attrs[source]):
        #    raise serializers.ValidationError('You can edit tags only for top level posts.')

        return attrs

    def validate_content(self, attrs, source):
        """
        Validation of `content` attribute.
        """
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

        if post.is_toplevel:
            if post.parent and not post.parent == post:
                raise serializers.ValidationError("The parent of a top level post must be either"
                                                  " emtpy or the post itself.")
        else:
            if not post.parent:
                raise serializers.ValidationError("Non top level posts must have parents.")

        if post.type == Post.ANSWER:
            if not post.parent.is_toplevel:
                raise serializers.ValidationError("Answers must have top level parents.")


#class ParentFieldValidator:
#    def __init__(self, serializer, attrs, source):
#        self.serializer = serializer
#        self.is_post_request = self.is_put_request = self.is_patch_request = False
#        #self.old_post = self.new_post = None
#        self.post = None
#
#        self._init_request_method()
#        #self._init_fields(attrs, source)
#
#    def _init_request_method(self):
#        method = self.context['request'].method.upper()
#        self.is_post_request = method == 'POST'
#        self.is_put_request = method == 'PUT'
#        self.is_patch_request = method == 'PATCH'
#
#    def _init_fields(self, attrs, source):
#        if self.is_patch_request or self.is_put_request:
#            self.post = self.serializer.object
#
#        if self.is_post_request:
#            self.post = Post()
#
#        for key, val in attrs.items():
#            setattr(self.post, key, val)
#
#        self.new_post.type = attrs.get('type', None) or self.new_post.type
#        self.new_post.is_toplevel = self.new_post.type in Post.TOP_LEVEL
#        self.new_post.parent = attrs.get('parent', None) or self.new_post.parent


## Custom permissions #############################################################################
class IsOwnerOrReadOnly(permissions.BasePermission):
    """
    Custom permission to allow only the author of a vote to delete it.
    """
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
    def has_object_permission(self, request, view, obj):
        # Read permissions are allowed to any request,
        # so we'll always allow GET, HEAD or OPTIONS requests.
        if request.method in permissions.SAFE_METHODS:
            return True

        # Write permissions are only allowed to the open posts.
        return obj.status == Post.OPEN