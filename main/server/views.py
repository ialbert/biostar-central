"""
Biostar views
"""
import difflib, time
from datetime import datetime, timedelta

from functools import partial
from collections import defaultdict
from main.server import html, models, const, formdef, action, notegen, auth
from main.server.html import get_page
from datetime import datetime

from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login, logout
from django.contrib import messages
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q
# the openid association model
from django_openid_auth.models import UserOpenID
from django.core.urlresolvers import reverse

# import all constants
from main.server.const import *

# activate logging
import logging
logger = logging.getLogger(__name__)


def update_counts(request, key, value):
    counts = request.session.get(SESSION_POST_COUNT,{})
    counts[key] = value
    request.session[SESSION_POST_COUNT] = counts

def get_post_manager(request):
    user = request.user
    if user.is_authenticated() and user.profile.can_moderate:
        return models.Post.objects
    else:
        return models.Post.open_posts
    
def index(request, tab="questions"):
    "Main page"
    
    user = request.user

    params = html.Params(tab=tab)
    
    # this will fill in the query (q) and the match (m)parameters
    params.parse(request)
    
    # update with counts
    counts = request.session.get(SESSION_POST_COUNT, {})

    # returns the object manager that contains all or only visible posts
    posts = get_post_manager(request)
    
    # apply search
    if params.q:
        pids  = action.search(params.q)
        posts = posts.filter(id__in=pids)
        tab   = '' # apply search sitewide rather than the selected tab
    
    # filter the posts by the tab that the user has selected
    if tab == "popular":
        posts = posts.filter(type=POST_QUESTION).order_by('-score')
    elif tab == "questions":
        posts = posts.filter(type=POST_QUESTION).order_by('-rank')
    elif tab == "unanswered":
        posts = posts.filter(type=POST_QUESTION, answer_count=0).order_by('-rank')
    elif tab == "recent":
        posts = posts.order_by('-creation_date')
    elif tab == 'planet':
        posts = posts.filter(type=POST_BLOG).order_by('-rank')
        models.decorate_posts(posts, user)
    elif tab == 'forum':
        posts = posts.filter(type=POST_FORUM).order_by('-rank')
    elif tab == 'tutorials':
        posts = posts.filter(type=POST_TUTORIAL).order_by('-rank')
    else:
        posts = posts.order_by('-rank')
    
    # reset the counts
    update_counts(request, tab, 0)

    page = get_page(request, posts, per_page=20)
    
    return html.template(request, name='index.html', page=page, params=params, counts=counts)

def show_tag(request, tag_name=None):
    "Display posts by a certain tag"
    params = html.Params(nav='', tab='')
    messages.warning(request, 'Filtering by tag: %s' % tag_name)
    posts = get_post_manager(request).filter(tag_set__name=tag_name)
    posts = posts.order_by('-creation_date')
    page  = get_page(request, posts, per_page=20)
    return html.template( request, name='index.html', page=page, params=params)

def show_user(request, uid, post_type=''):
    "Displays posts by a user"

    user = models.User.objects.filter(id=uid).select_related('profile').all()[0]
    params = html.Params(nav='', tab='')

    # notification
    messages.info(request, 'Filtering by user: %s' % user.profile.display_name)
   
    post_type = POST_REV_MAP.get(post_type.lower())
    if post_type:
        posts = get_post_manager(request).filter(type=post_type, author=user).order_by('-creation_date')
    else:
        posts = get_post_manager(request).filter(type__in=POST_TOPLEVEL, author=user).order_by('-creation_date')
    page  = get_page(request, posts, per_page=20)
    return html.template( request, name='index.html', page=page, params=params)


def user_profile(request, uid, tab='activity'):
    "User's profile page"

    user = request.user
    target = models.User.objects.get(id=uid)
    awards = []
    page   = None
    
    # some information is only visible to the user
    target.writeable = auth.authorize_user_edit(target=target, user=user, strict=False)
    target.showall = (target == user)

    params = html.Params(tab=tab)

    # these do not actually get executed unless explicitly rendered in the page
    bookmarks = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).select_related('post', 'post__author__profile').order_by('id')
    awards = models.Award.objects.filter(user=target).select_related('badge').order_by('-date')

 # we need to collate and count the awards
    answer_count = models.Post.objects.filter(author=target, type=POST_ANSWER).count()
    question_count = models.Post.objects.filter(author=target, type=POST_QUESTION).count()
    comment_count = models.Post.objects.filter(author=target, type=POST_COMMENT).count()
    post_count = models.Post.objects.filter(author=target).count()
    vote_count = models.Vote.objects.filter(author=target).count()
    award_count = models.Award.objects.filter(user=target).count()
    note_count  = models.Note.objects.filter(target=target, unread=True).count()
    bookmarks_count  = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).count()
    
    
    if tab == 'activity':
        notes = models.Note.objects.filter(target=target).select_related('author', 'author__profile', 'root').order_by('-date')
        page  = get_page(request, notes, per_page=15)
        # we evalute it here so that subsequent status updates won't interfere
        page.object_list = list(page.object_list)
        if user==target:
            models.Note.objects.filter(target=target).update(unread=False)
            models.UserProfile.objects.filter(user=target).update(new_messages=0)
            note_count = 0
        
    elif tab == 'bookmarks':
        bookmarks = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).select_related('post', 'post__author__profile').order_by('id')
        page  = get_page(request, bookmarks, per_page=5)
    
    elif tab =="moderator":
       
        notes = models.Note.objects.filter(target=target, type=NOTE_MODERATOR).select_related('author', 'author__profile', 'root').order_by('-date')
        page  = get_page(request, notes, per_page=15)
        print 'HERE'

    params.update(dict(question_count=question_count, answer_count=answer_count, note_count=note_count, bookmarks_count=bookmarks_count,
            comment_count=comment_count, post_count=post_count, vote_count=vote_count, award_count=award_count))
    
    return html.template(request, name='user.profile.html', awards=awards,
        user=request.user,target=target, params=params, page=page)

def user_list(request):
    search  = request.GET.get('m','')[:80] # trim for sanity
    params = html.Params(nav='users')
    if search:
        query = Q(profile__display_name__icontains=search)
        users = models.User.objects.filter(query).select_related('profile').order_by("-profile__score")
    else:
        users = models.User.objects.select_related('profile').order_by("-profile__score")
    page  = get_page(request, users, per_page=24)
    return html.template(request, name='user.list.html', page=page, params=params)

def tag_list(request):
    tags = models.Tag.objects.all().order_by('-count')
    page = get_page(request, tags, per_page=50)
    params = html.Params(nav='tags')
    return html.template(request, name='tag.list.html', page=page, params=params)

def badge_list(request):
    badges = models.Badge.objects.filter(secret=False).order_by('-count', '-type')
    params = html.Params(nav='badges')
    return html.template(request, name='badge.list.html', badges=badges, params=params)
 
def post_show(request, pid):
    "Returns a question with all answers"
    user = request.user

    query = get_post_manager(request)

    try:
        root = query.get(id=pid)
        # update the views for the question
        root.update_views(request)
    except models.Post.DoesNotExist, exc:
        messages.warning(request, 'The post that you are looking for does not exists. Perhaps it was deleted!')
        return html.redirect("/")
    
    # get all answers to the root
    children = models.Post.objects.filter(root=root).exclude(type=POST_COMMENT).select_related('author', 'author__profile').order_by('-accepted', '-score')
    
    # comments need to be displayed by creation date
    comments = models.Post.objects.filter(root=root, type=POST_COMMENT).select_related('author', 'author__profile').order_by('creation_date')

    all = [ root ] + list(children) + list(comments)
    # add the various decorators
    models.decorate_posts(all, user)
    
    # these are all the answers
    answers = [ o for o in children if o.type == POST_ANSWER ]
    
    # get all the comments
    tree = defaultdict(list)
    for comment in comments:  
        tree[comment.parent_id].append(comment)
   
    # generate the tag cloud
    #tags = models.Tag.objects.all().order_by('-count')[:50]
    
    return html.template( request, name='post.show.html', root=root, answers=answers, tree=tree )
 
def post_redirect(post):
    return html.redirect( post.get_absolute_url() )
    
@login_required(redirect_field_name='/openid/login/')
def new_comment(request, pid=0):
    "Shortcut to new comments"
    return new_post(request=request, pid=pid, post_type=POST_COMMENT)

@login_required(redirect_field_name='/openid/login/')
def new_answer(request, pid):
    return new_post(request=request, pid=pid, post_type=POST_ANSWER)
    
@login_required(redirect_field_name='/openid/login/')
def new_post(request, pid=0, post_type=POST_QUESTION):
    "Handles the creation of a new post"
    
    user   = request.user
    name   = "post.edit.html"
    parent = models.Post.objects.get(pk=pid) if pid else None
    root   = parent.root if parent else None
    toplevel = (pid == 0)
    factory  = formdef.ChildContent if pid else formdef.TopLevelContent
    
    params = html.Params(tab='new', title="New post", toplevel=toplevel)
    
    if request.method == 'GET':
        # no incoming data, render form
        form = factory()
        return html.template(request, name=name, form=form, params=params)
    
    # process the incoming data
    assert request.method == 'POST', "Method=%s" % request.method
    
    form = factory(request.POST)
    if not form.is_valid():
        # returns with an error message
        return html.template(request, name=name, form=form, params=params)

    # form is valid at this point, create the post
    params = dict(author=user, type=post_type, parent=parent, root=root)
    params.update(form.cleaned_data)
    
    with transaction.commit_on_success():
        post = models.Post.objects.create(**params)
        post.set_tags()
        post.save()
    return post_redirect(post)

@login_required(redirect_field_name='/openid/login/')
def post_edit(request, pid=0):
    "Handles the editing of an existing post"
    
    user = request.user
    name   = "post.edit.html"
    post = models.Post.objects.get(pk=pid)

    if not post.open and not user.can_moderate:
        messages.error(request, 'Post is closed. It may not be edited.')
        return post_redirect(post.root)
    
    # verify that this user may indeed modify the post
    auth.authorize_post_edit(post=post, user=request.user, strict=True)
    
    toplevel = post.top_level
    factory  = formdef.TopLevelContent if toplevel else formdef.ChildContent

    params = html.Params(tab='edit', title="Edit post", toplevel=toplevel)
    if request.method == 'GET':
        # no incoming data, render prefilled form
        form = factory(initial=dict(title=post.title, content=post.content, tag_val=post.tag_val, type=post.type))
        return html.template(request, name=name, form=form, params=params)

    # process the incoming data
    assert request.method == 'POST', "Method=%s" % request.method
    form = factory(request.POST)
    if not form.is_valid():
        # returns with an error message
        return html.template(request, name=name, form=form, params=params)

    # form is valid now set the attributes
    for key, value in form.cleaned_data.items():
        setattr(post, key, value)
        post.set_tags()
    models.create_revision(post)
    return post_redirect(post)
    
def revision_show(request, pid):
    post = models.Post.objects.get(pk=pid)
    revs = post.revisions.order_by('-date').select_related('author')
    return html.template(request, name='revision.show.html', revs=revs, post=post)
   
@login_required(redirect_field_name='/openid/login/')
def add_comment(request, pid):
    "Adds a comment"
    parent  = models.Post.objects.get(pk=pid)
    content = request.POST['text'].strip()
    if len(content)<1:
        messages.warning(request, 'Comment too short!')
        return post_redirect(parent)
        
    comment = models.Post(author=request.user, parent=parent, post_type=POST_COMMENT, creation_date=datetime.now())
    comment.save()
    comment.create_revision(content=content)
    
    return post_redirect(comment)

def blog_redirect(request, pid):
    "Used to be able to count the views for a blog"
    blog = models.Post.objects.get(id=pid, type=POST_BLOG)
    if blog.update_views(request):
        blog.rank_change(sign=1, hours=6) 
        blog.save()
    return html.redirect( blog.get_absolute_url() )

def modlog_list(request):
    "Lists of all moderator actions"
    mods = models.Note.objects.filter(type=NOTE_MODERATOR).select_related('sender', 'target', 'post', 'sender_profile').order_by('-date')
    page = get_page(request, mods)
    return html.template(request, name='modlog.list.html', page=page)

