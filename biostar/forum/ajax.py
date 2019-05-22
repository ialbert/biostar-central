
from django.contrib import messages
from django.contrib.auth.decorators import login_required

from django.shortcuts import render

from biostar.forum import forms, auth
from biostar.forum.models import Post, Vote
from biostar.utils.decorators import ajax_error, ajax_error_wrapper, ajax_success


@ajax_error_wrapper(method="POST")
def ajax_vote(request):
    user = request.user
    type_map = dict(upvote=Vote.UP, bookmark=Vote.BOOKMARK, accept=Vote.ACCEPT)

    vote_type = request.POST.get('vote_type')
    vote_type = type_map.get(vote_type)
    post_uid = request.POST.get('post_uid')

    # Check the post that is voted on.
    post = Post.objects.filter(uid=post_uid).first()

    if post.author == user and vote_type == Vote.UP:
        return ajax_error("You can not upvote your own post.")

    if post.author == user and vote_type == Vote.ACCEPT:
        return ajax_error("You can not accept your own post.")

    if post.root.author != user and vote_type == Vote.ACCEPT:
        return ajax_error("Only the person asking the question may accept this answer.")

    msg, vote = auth.apply_vote(post=post, user=user, vote_type=vote_type)
    return ajax_success(msg=msg)


def ajax_test(request):
    """
    Creates a commment on a top level post.
    """
    msg="OK"
    print (f"HeRe= {request.POST} ")
    return ajax_success(msg=msg)
