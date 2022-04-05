
from pagedown.widgets import PagedownWidget
import os
import langdetect
import logging
from functools import wraps
from django import forms
from django.utils.safestring import mark_safe
from django.core.exceptions import ValidationError
from django.shortcuts import reverse
from django.template import loader
from django.contrib import messages
from django.shortcuts import redirect, render
from django.contrib.auth.decorators import login_required
from django.conf import settings
from biostar.accounts.views import user_moderate as account_moderate
from biostar.accounts.models import Profile, User
from biostar.utils.decorators import check_params
from biostar.forum.models import Post, delete_post_cache, Log
from biostar.forum import auth, const, util


logger = logging.getLogger('engine')

def post_exists(func):
    """
    Ensure uid passed to view function exists.
    """

    @wraps(func)
    def _wrapper_(request, **kwargs):
        uid = kwargs.get('uid')
        post = Post.objects.filter(uid=uid).exists()
        if not post:
            messages.error(request, "Post does not exist.")
            return redirect(reverse("post_list"))
        return func(request, **kwargs)

    return _wrapper_


class PostModForm(forms.Form):

    def __init__(self, post, request, user, *args, **kwargs):
        self.post = post
        self.user = user
        self.request = request

        super(PostModForm, self).__init__(*args, **kwargs)

        choices = [
            ('delete', "Delete post"),
            ('open', "Open post"),
        ]

        # Top level posts may be bumped.
        if post.is_toplevel:
            choices += [('bump', "Bump post")]
        elif post.is_comment or post.is_answer:
            choices += [('relocate', "Move post")]

        # Punitive options.
        choices.extend((
                ('offtopic', "Mark as offtopic"),
                ('spam', "Mark as spam (suspend author)"),

        ))

        if settings.ALLOW_POST_CLOSING:
            choices.extend(('close', "Close post"))

        self.fields['action'] = forms.CharField(widget=forms.RadioSelect(choices=choices), required=True)

    def clean(self):

        if not self.user.profile.is_moderator:
            raise forms.ValidationError("You need to be a moderator to preform that action.")

        return self.cleaned_data


@check_params(allowed=const.ALLOWED_PARAMS)
@login_required
def user_moderate(request, uid):

    def callback():
        source = request.user
        target = User.objects.filter(id=uid).first()
        text = f"changed user state to {target.profile.get_state_display()}"
        auth.db_logger(user=source, text=text, target=target)
        return

    result = account_moderate(request, uid=uid, callback=callback)

    return result


@check_params(allowed=const.ALLOWED_PARAMS)
@post_exists
@login_required
def post_moderate(request, uid):
    """Used to make display post moderate form given a post request."""

    user = request.user
    post = Post.objects.filter(uid=uid).first()

    if request.method == "POST":
        form = PostModForm(post=post, data=request.POST, user=user, request=request)
        if form.is_valid():
            action = form.cleaned_data.get('action')
            url = moderate(request=request, post=post, action=action)
            return redirect(url)
        else:
            errors = ','.join([err for err in form.non_field_errors()])
            messages.error(request, errors)
            return redirect(reverse("post_view", kwargs=dict(uid=post.root.uid)))
    else:
        form = PostModForm(post=post, user=user, request=request)

    context = dict(form=form, post=post, user=user, ALLOW_POST_CLOSING=settings.ALLOW_POST_CLOSING)
    return render(request, "forms/form_moderate.html", context)


def mod_rationale(post, user, template, ptype=Post.ANSWER, extra_context=dict()):
    tmpl = loader.get_template(template)
    context = dict(user=post.author)
    context.update(extra_context)
    content = tmpl.render(context)

    # Load answer explaining post being off topic.
    post = Post.objects.create(content=content, type=ptype, parent=post, root=post.root, author=user)

    return post


def removal_condition(post, user, age=1):
    """
    Removal condition for the post.
    """

    # Only authors may remove their own posts
    if post.author != user:
        return False

    # Post older than a day may not be removed
    if post.age_in_days > age:
        return False

    # If the post has children it may not be removed
    if Post.objects.filter(parent=post).exclude(pk=post.id):
        return False

    # If the post has votes it may not be removed
    if post.vote_count:
        return False

    return True


def delete_post(request, post, **kwargs):
    """
    Post may be marked as deleted or removed entirely
    """
    user = request.user

    # Decide on the removal
    remove = removal_condition(post, user)

    if remove:
        msg = f"removed post"
        messages.info(request, mark_safe(msg))
        auth.db_logger(user=user, post=post, text=msg)
        post.delete()
        # Deleted children should return root url.
        url = "/" if post.is_toplevel else post.root.get_absolute_url()
    else:
        Post.objects.filter(uid=post.uid).update(status=Post.DELETED)
        post.recompute_scores()
        msg = f"deleted post"
        messages.info(request, mark_safe(msg))
        auth.db_logger(user=user, post=post, text=msg)
        url = post.get_absolute_url()

    # Recompute post score.
    if not post.is_toplevel:
        post.root.recompute_scores()

    return url


def open(request, post, **kwargs):
    if post.is_spam and post.author.profile.low_rep:
        post.author.profile.bump_over_threshold()

    user = request.user
    Post.objects.filter(uid=post.uid).update(status=Post.OPEN, spam=Post.NOT_SPAM)
    post.recompute_scores()

    post.root.recompute_scores()
    msg = f"opened post"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    auth.db_logger(user=user, text=f"{msg}", post=post)
    return url


def bump(request, post, **kwargs):
    now = util.now()
    user = request.user

    Post.objects.filter(uid=post.uid).update(lastedit_date=now, rank=now.timestamp())
    msg = f"bumped post"
    url = post.get_absolute_url()
    messages.info(request, mark_safe(msg))
    auth.db_logger(user=user, text=f"{msg}", post=post)

    return url


def change_user_state(mod, target, state):
    """
    Changes user state.
    """

    # Only moderators may change user states.
    if not mod.profile.is_moderator:
        logger.error(f"{mod} is not a moderator")
        return

    # Cannot moderate self.
    if mod == target:
        logger.error(f"{mod} cannot moderate themselves")
        return

    # The target may not be a moderator.
    if target.profile.is_moderator:
        logger.info(f"{mod} cannot alter state on a moderator {target}")
        return

    # Set the new state.
    target.profile.state = state
    target.profile.save()

    # Generate the logging message.
    msg = f"changed user state to {target.profile.get_state_display()}"
    auth.db_logger(user=mod, action=Log.MODERATE, target=target, text=msg, post=None)


def toggle_spam(request, post, **kwargs):
    """
    Toggles spam status on post based on a status
    """

    url = post.get_absolute_url()

    # Moderators may only be suspended by admins (TODO).
    if post.author.profile.is_moderator and post.spam in (Post.DEFAULT, Post.NOT_SPAM):
        messages.warning(request, "cannot toggle spam on a post created by a moderator")
        return url

    # The user performing the action.
    user = request.user

    # Drop the cache for the post.
    delete_post_cache(post)

    # Current state of the toggle.
    if post.is_spam:
        Post.objects.filter(id=post.id).update(spam=Post.NOT_SPAM, status=Post.OPEN)
    else:
        Post.objects.filter(id=post.id).update(spam=Post.SPAM, status=Post.CLOSED)

    # Refetch up to date state of the post.
    post = Post.objects.filter(id=post.id).get()

    # Set the state for the user (only non moderators are affected)
    state = Profile.SUSPENDED if post.is_spam else Profile.NEW

    # Apply the user change.
    change_user_state(mod=user, target=post.author, state=state)

    # Generate logging messages.
    if post.is_spam:
        text = f"marked post as spam"
    else:
        text = f"restored post from spam"

        # Set indexed flag to False, so it's removed from spam index
        Post.objects.filter(id=post.id).update(indexed=False)

    # Set a logging message.
    messages.success(request, text)

    # Submit the log into the database.
    auth.db_logger(user=user, action=Log.MODERATE, target=post.author, text=text, post=post)

    url = post.get_absolute_url()

    return url


def close(request, post, **kwargs):

    """
    Close this post and provide a rationale for closing as well.
    """
    user = request.user
    Post.objects.filter(uid=post.uid).update(status=Post.CLOSED)
    # Generate a rationale post on why this post is closed.
    rationale = mod_rationale(post=post, user=user,
                              template="messages/closed.md")
    msg = "closed"
    url = rationale.get_absolute_url()
    messages.info(request, mark_safe(msg))
    auth.db_logger(user=user, text=f"{msg}", post=post)
    return url


def off_topic(request, post, **kwargs):
    """
    Marks post as off topic. Generate off topic comment.
    """
    user = request.user
    if "offtopic" not in post.tag_val:
        post.tag_val += ",offtopic "
        post.save()

        # Generate off comment.
        template = 'messages/offtopic.md'
        tmpl = loader.get_template(template)
        context = dict(post=post)
        content = tmpl.render(context)

        auth.create_post(ptype=Post.COMMENT, parent=post, content=content, title='', author=request.user, request=request)
        msg = "off topic"
        messages.info(request, mark_safe(msg))
        auth.db_logger(user=user, text=f"{msg}", post=post)
    else:
        messages.warning(request, "post has been already tagged as off topic")

    url = post.get_absolute_url()
    return url


def relocate(request, post, **kwds):
    """
    Moves an answer to a comment and viceversa.
    """
    url = post.get_absolute_url()

    if post.is_toplevel:
        messages.warning(request, "cannot relocate a top level post")
        return url

    if post.type == Post.COMMENT:
        msg = f"relocated comment to answer"
        post.type = Post.ANSWER
    else:
        msg = f"relocated answer to comment"
        post.type = Post.COMMENT

    post.parent = post.root
    post.save()
    post.update_parent_counts()

    auth.db_logger(user=request.user, post=post, text=f"{msg}")
    messages.info(request, msg)
    return url


def moderate(request, post, action):

    # Bind an action to a function.
    action_map = {
        "spam": toggle_spam,
        "bump": bump,
        "open": open,
        "offtopic": off_topic,
        "delete": delete_post,
        "close": close,
        "relocate": relocate,
    }

    if action in action_map:
        mod_func = action_map[action]
        url = mod_func(request=request, post=post)
    else:
        url = post.get_absolute_url()
        msg = "Unknown moderation action given."
        logger.error(msg)

    return url

