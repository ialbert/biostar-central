__author__ = 'ialbert'
import json, traceback, logging
from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User
from django.http import HttpResponse
from functools import partial
from django.db import transaction
from django.db.models import F


def json_response(adict, **kwd):
    """Returns a http response in JSON format from a dictionary"""
    return HttpResponse(json.dumps(adict), **kwd)


logger = logging.getLogger(__name__)


def ajax_msg(msg, status, **kwargs):
    payload = dict(status=status, msg=msg)
    payload.update(kwargs)
    return json_response(payload)


ajax_success = partial(ajax_msg, status='success')
ajax_error = partial(ajax_msg, status='error')


class ajax_error_wrapper(object):
    """
    Used as decorator to trap/display  errors in the ajax calls
    """

    def __init__(self, f):
        self.f = f

    def __call__(self, request):
        try:
            if request.method != 'POST':
                return ajax_error('POST method must be used.')

            if not request.user.is_authenticated():
                return ajax_error('You must be logged in to do that')

            value = self.f(request)
            return value
        except Exception, exc:
            traceback.print_exc()
            return ajax_error('Error: %s' % exc)


POST_TYPE_MAP = dict(vote=Vote.UP, bookmark=Vote.BOOKMARK, accept=Vote.ACCEPT)


@transaction.atomic
def perform_vote(post, user, vote_type):
    """
    Stores a vote and run side effects like updating the author reputation.
    """
    # Only maintain one vote for each user/post pair.
    votes = Vote.objects.filter(author=user, post=post, type=vote_type)
    if votes:
        vote = votes[0]
        msg = "%s removed" % vote.get_type_display()
        change = -1
    else:
        change = +1
        vote = Vote.objects.create(author=user, post=post, type=vote_type)
        msg = "%s added" % vote.get_type_display()

    if post.author != user:
        # Update the user reputation only if the author is different.
        User.objects.filter(pk=post.author.id).update(score=F('score') + change)

    # The thread score represents all votes in a thread
    Post.objects.filter(pk=post.root_id).update(thread_score=F('thread_score') + change)

    if vote.type == Vote.BOOKMARK:
        # Apply the vote
        Post.objects.filter(pk=post.id).update(book_count=F('book_count') + change, vote_count=F('vote_count') + change)
        Post.objects.filter(pk=post.id).update(subs_count=F('subs_count') + change)
        Post.objects.filter(pk=post.root_id).update(subs_count=F('subs_count') + change)

    elif vote_type == Vote.ACCEPT:
        if change > 0:
            # There does not seem to be a negation operator for F objects.
            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=True)
            Post.objects.filter(pk=post.root_id).update(has_accepted=True)
        else:
            Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change, has_accepted=False)
            Post.objects.filter(pk=post.root_id).update(has_accepted=False)
    else:
        Post.objects.filter(pk=post.id).update(vote_count=F('vote_count') + change)

    # Clear old votes.
    if votes:
        votes.delete()

    return msg


@ajax_error_wrapper
def vote_handler(request):
    """
    Handles all voting on posts via AJAX.
    """
    user = request.user
    vote_type = request.POST['vote_type']
    vote_type = POST_TYPE_MAP[vote_type]
    post_id = request.POST['post_id']
    post = Post.objects.get(pk=post_id)

    # Create a new vote for validation (not storing it to the db yet).
    vote = Vote(author=user, post=post, type=vote_type)

    # Perform the vote validation.
    validation_result = validate_vote(vote)

    # Parse the result fo the validation and return a proper ajax message.
    if validation_result == VALID_VOTE:
        with transaction.atomic():
            msg = perform_vote(post=post, user=user, vote_type=vote_type)
        return ajax_success(msg)

    if validation_result == UPVOTED_OWN_POST:
        return ajax_error("You can't upvote your own post.")

    if validation_result == ACCEPTED_NOT_OWN_QUESTION:
        return ajax_error("Only the person asking the question may accept this answer.")

    return ajax_error("{}.".format(VOTE_VALIDATION_MSGS[validation_result]))


VALID_VOTE = 0
DOWNVOTE = 1
UPVOTED_OWN_POST = 2
ACCEPTED_NOT_ANSWER = 3
ACCEPTED_NOT_OWN_QUESTION = 4
VOTE_VALIDATION_MSGS = {
    VALID_VOTE: 'Ok',
    DOWNVOTE: 'Downvotes are not allowed',
    UPVOTED_OWN_POST: 'You are not allowed to upvote your own posts',
    ACCEPTED_NOT_ANSWER: 'Only answer posts can be accepted',
    ACCEPTED_NOT_OWN_QUESTION: 'Only the author of a question can accept a relative answer',
}


def validate_vote(vote):
    """
    Vote validation based on some rules.
    """
    # Rule 1: downvotes not allowed.
    if vote.type == Vote.DOWN:
        return DOWNVOTE

    # Rule 2: a user can not upvote her own post.
    if (vote.type == Vote.UP and
        vote.author == vote.post.author):
        return UPVOTED_OWN_POST

    # Rule 3: accept votes are only allowed for posts of type "answer".
    if (vote.type == Vote.ACCEPT and
        not vote.post.type == Post.ANSWER):
        return ACCEPTED_NOT_ANSWER

    # Rule 4: the author of an accept vote must match the author of the root post (the
    # original question).
    if (vote.type == Vote.ACCEPT and
        not vote.author == vote.post.root.author):
        return ACCEPTED_NOT_OWN_QUESTION

    # Rule 5: users can not accept their own answers.
    #if (vote.type == Vote.ACCEPT and
    #    vote.author == vote.post.author):
    #    return ...

    return VALID_VOTE