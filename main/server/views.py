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

def get_posts(request):
    "Returns a common queryset that can be used to select questions"
    if request.user.can_moderate:
        query = models.Post.all_posts
    else:
        query = models.Post.open_posts
    return query

def update_counts(request, key, value):
    counts = request.session.get(SESSION_POST_COUNT,{})
    counts[key] = value
    request.session[SESSION_POST_COUNT] = counts

def index(request, tab="questions"):
    "Main page"
    
    params = html.Params(tab=tab)
    
    # this will fill in the query (q) and the match (m)parameters
    params.parse(request)
    
    # update with counts
    counts = request.session.get(SESSION_POST_COUNT, {})

    # returns the object manager that contains all or only visible posts
    posts = get_posts(request)
    
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
    elif tab == 'forum':
        posts = posts.filter(type=POST_FORUM).order_by('-rank')
    elif tab == 'guides':
        posts = posts.filter(type=POST_GUIDE).order_by('-rank')
    else:
        posts = posts.order_by('-rank')
    
    # reset the counts
    update_counts(request, tab, 0)

    page = get_page(request, posts, per_page=20)
    
    return html.template(request, name='index.html', page=page, params=params, counts=counts)

def show_tag(request, tag_name=None):
    "Display posts by a certain tag"
    params = html.Params(nav='', tab='')
    params.setr('Filtering by tag: %s' % tag_name)
    posts = get_posts(request).filter(tag_set__name=tag_name)
    posts = posts.order_by('-creation_date')
    page  = get_page(request, posts, per_page=20)
    return html.template( request, name='index.html', page=page, params=params)

def show_user(request, uid, post_type=''):
    "Displays posts by a user"

    user = models.User.objects.filter(id=uid).select_related('profile').all()[0]
    params = html.Params(nav='', tab='')
    params.setr('Filtering by user: %s' % user.profile.display_name)
    post_type = POST_REV_MAP.get(post_type.lower())
    if post_type:
        posts = get_posts(request).filter(type=post_type, author=user).order_by('-creation_date')
    else:
        posts = get_posts(request).filter(type__in=POST_TOPLEVEL, author=user).order_by('-creation_date')
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

def question_unanswered(request, uid=0, post_type=None):
    "Lists all the questions"
    params = html.Params()
    params.setr('Filter: unanswered')
    qs = get_posts(request).filter(answer_count=0, type=POST_QUESTION)
    page = get_page(request, qs) 
    return html.template(request, name='post.list.html', page=page, params=params)

def question_tagged(request, tag_name):
    params = html.Params()
    params.setr('Tag: %s' % tag_name)
    qs = get_posts(request).filter(tag_set__name=tag_name)
    page = get_page(request, qs) 
    return html.template(request, name='post.list.html', page=page, params=params)
  
def post_show(request, pid):
    "Returns a question with all answers"
    user = request.user

    query = models.Post.objects
    try:
        root = query.get(id=pid)
        # update the views for the question
        root.update_views(request)
        auth.authorize_post_edit(post=root, user=request.user, strict=False)
    except models.Post.DoesNotExist, exc:
        messages.warning(request, 'The post that you are looking for does not exists. Perhaps it was deleted!')
        return html.redirect("/")
    
    # get all answers to the root
    children = models.Post.objects.filter(root=root).select_related('author', 'author__profile').order_by('-accepted', '-score')
   
    # these are all the answers
    answers = [ o for o in children if o.type == POST_ANSWER ]
   
    # all objects with votes
    all = list(children) + [ root ]
        
    if request.user.is_authenticated():
        votes = models.Vote.objects.filter(author=request.user, post__id__in = [ p.id for p in all ] ) 
        up_votes  = set(vote.post.id for vote in votes if vote.type == const.VOTE_UP)
        bookmarks = set(vote.post.id for vote in votes if vote.type == const.VOTE_BOOKMARK)
    else:
        up_votes = down_votes = bookmarks = set()

    # decorate the posts with extra attributes for easier rendering
    for post in all :
        post.writeable = auth.authorize_post_edit(post=post, user=request.user, strict=False)
        post.upvoted   = post.id in up_votes
        post.bookmarked = post.id in bookmarks
    
    # get all the comments
    tree = defaultdict(list)
    comments = models.Post.objects.filter(root=root, type=POST_COMMENT).select_related('author', 'author__profile').order_by('creation_date')
    for comment in comments:
        comment.writeable = auth.authorize_post_edit(post=comment, user=request.user, strict=False)
        comment.upvoted   = comment.id in up_votes
        tree[comment.parent_id].append(comment)
   
    # generate the tag cloud
    tags = models.Tag.objects.all().order_by('-count')[:50]
    
    return html.template( request, name='post.show.html', root=root, answers=answers, tree=tree, tags=tags )
 
def post_redirect(post, anchor=None):
    """
    Shows a post in full context
    """
    # get the root of a post
    root = post.root or post
    pid, slug = root.id, root.slug
    anchor = anchor or post.id
    url = '/post/show/%s/%s/#%s' % (pid, slug, anchor)
    return html.redirect(url)
    
@login_required(redirect_field_name='/openid/login/')
def new_comment(request, parentid=0):
    "Shortcut to new comments"
    return post_edit(request=request, pid=0, parentid=parentid, post_type=POST_COMMENT)

@login_required(redirect_field_name='/openid/login/')
def new_answer(request,  parentid=0):
    "Shortcut to new answers"
    return post_edit(request=request, pid=0, parentid=parentid, post_type=POST_ANSWER)

@login_required(redirect_field_name='/openid/login/')
def new_question(request, pid=0):
    "Shortcut to new questions"
    return post_edit(request=request, pid=0, parentid=0, post_type=POST_QUESTION)

@login_required(redirect_field_name='/openid/login/')
def post_edit(request, pid=0, parentid=0, post_type=POST_QUESTION):
    """
    Handles post related edits for all posts
    """
    
    # a shortcut
    user = request.user

    # new post creation
    newpost = (pid == 0)

    # incoming data in the request 
    form_data = (request.method == 'POST')
    
    # sanity check
    assert post_type in POST_MAP, 'Invalid post_type %s' % post_type
  
    # select the form factory from the post types
    use_post_form = (post_type in POST_TOPLEVEL)

    if use_post_form:
        factory = formdef.PostForm
    else:
        factory = formdef.ContentForm
    
    # find the parent if it exists
    if parentid:
        parent = models.Post.objects.get(pk=parentid)
        root   = parent.root or parent
    else:
        parent = root = None

    # this is the template name
    tmpl_name = "post.edit.html"
    
    # deal with new post creation first
    if newpost:
        # this here is to customize the output
        params = html.Params(tab='new', title="New post", use_post_form=use_post_form ) 

        if form_data:
            form = factory(request.POST)
            if not form.is_valid():
                return html.template(request, name=tmpl_name, form=form, params=params)
            params = dict(author=user, type=post_type, parent=parent, root=root, creation_date=datetime.now())
            params.update(form.cleaned_data)            
            with transaction.commit_on_success():
                post = models.Post.objects.create(**params)
                post.set_tags()
                models.create_revision(post)
            return post_redirect(post)
        else:
            form = factory()
            return html.template( request, name=tmpl_name, form=form, params=params)

    #
    # at this point we are dealing with a post editing action
    #
    assert pid, 'Only post modification should follow after this point'
    
    # when we edit a post we keep the original post type (for now)
    post = models.Post.objects.get(pk=pid)
    parent = post.parent
    post_type = post.type
   
    # select the form factory from the post types
    use_post_form = (post_type in POST_TOPLEVEL)
    if use_post_form:
        factory = formdef.PostForm
    else:
        factory = formdef.ContentForm

    # verify that this user may indeed modify the post
    auth.authorize_post_edit(post=post, user=request.user, strict=True)
    
    params = html.Params(title="Edit %s" % post.get_type_display(), use_post_form=use_post_form, tab=None) 

    # no form data coming, return the editing form
    if not form_data:
        form = factory(initial=dict(title=post.title, content=post.content, tag_val=post.tag_val))
        return html.template(request, name=tmpl_name, form=form, params=params)
    else:
        # we have incoming form data for posts
        form = factory(request.POST)
        if not form.is_valid():
            return html.template(request, name=tmpl_name, form=form, params=params)
         
        with transaction.commit_on_success():
            for key, value in form.cleaned_data.items():
                setattr(post, key, value)
            post.set_tags()
            models.create_revision(post)

        return post_redirect(post)    

def revision_list(request, pid):
    post = models.Post.objects.get(pk=pid)
    revs = post.revisions.order_by('-date').select_related('author')
    return html.template(request, name='revision.list.html', revs=revs, post=post)
   
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

#
# helper methods for json returns 
#
def ajax_msg(msg, status):
    return html.json_response(dict(status=status, msg=msg))
    
ajax_success = partial(ajax_msg, status='success')
ajax_error   = partial(ajax_msg, status='error')

class ajax_error_wrapper(object):
    "used as decorator to trap/display  errors in the ajax calls"
    def __init__(self, f):
        self.f = f
        
    def __call__(self, *args, **kwds):
        try:
            value = self.f(*args, **kwds)
            return value
        except Exception,exc:
            return ajax_error('Error: %s' % exc)

@ajax_error_wrapper           
def vote(request):
    "Handles all voting on posts"
    
    if request.method != 'POST':
        return ajax_error('POST method must be used')
        
    author = request.user
    if not author.is_authenticated():
        return ajax_error('You must be logged in to vote')
            
    # attempt to find the post and vote
    post_id = int(request.POST.get('post'))
    post = models.Post.objects.get(id=post_id)

    # get vote type
    type = request.POST.get('type')
    
    # remap to actual type
    type = dict(upvote=VOTE_UP, accept=VOTE_ACCEPT, bookmark=VOTE_BOOKMARK).get(type)
        
    if not type:
        return ajax_error('invalid vote type')
            
    if type == VOTE_UP and post.author == author:
        return ajax_error('You may not vote on your own post')
    
    if type == VOTE_ACCEPT and post.root.author != author:
        return ajax_error('Only the original poster may accept an answer')
        
    # see if there is an existing vote of this type
    old_vote = post.get_vote(author, type)

    if old_vote:
        msg = '%s removed' % old_vote.get_type_display()
        old_vote.delete()
        logger.info('%s\t%s\t%s' % (author.id, post.id, msg) )
        return ajax_success(msg)
    
    if type == VOTE_BOOKMARK:
        vote = post.add_vote(author, type)
        return ajax_success('%s added' % vote.get_type_display())
    
    today = datetime.now()
    shift = timedelta(seconds=VOTE_SESSION_LENGTH)
    past  = today - shift
    count = models.Vote.objects.filter(author=author, date__gt=past).count()
    avail  = MAX_VOTES_PER_SESSION - count
    
    if avail <= 0:
        msg = "You ran out of votes ;-) there will more in a little while"
        logger.info('%s\t%s\t%s' % (author.id, post.id, "out of votes") )
        return ajax_error(msg)
    else:
        # log all voting into the server log
        vote = post.add_vote(author, type)
        msg  = '%s added. %d votes left' % (vote.get_type_display(), avail-1)
        logger.info('%s\t%s\t%s' % (author.id, post.id, msg) )
        return ajax_success(msg)

@login_required(redirect_field_name='/openid/login/')
def moderate_post(request, pid, action):
    
    user = request.user
    if request.method != 'POST':
        return html.json_response({'status':'error', 'msg':'Only POST requests are allowed'})        
   
    post = models.Post.objects.get(id=pid)
    status = auth.authorize_post_edit(user=user, post=post, strict=False)
    if not status:
        return html.json_response({'status':'error', 'msg':'You do not have permission to moderate this post.'})        
        
    action_map = { 'close':models.REV_CLOSE, 'reopen':models.REV_REOPEN,
                  'delete':models.REV_DELETE, 'undelete':models.REV_UNDELETE }
    
    models.moderator_action(post=post, action=action_map[action], user=user)
    return html.json_response({'status':'success', 'msg':'%s performed' % action})
        
@login_required(redirect_field_name='/openid/login/')
def moderate_user(request, uid, action):

    if request.method != 'POST':
        return html.json_response({'status':'error', 'msg':'Only POST requests are allowed'})        
    
    moderator = request.user
    target = models.User.objects.get(id=uid)

    if not target.profile.authorize(moderator=moderator):
        return html.json_response({'status':'error', 'msg':'You do not have permission to moderate this user.'})        
    
    if action == 'suspend':
        target.profile.suspended = True
        target.profile.save()
        text = notegen.suspend(target)
        models.Note.send(target=moderator, content=text, sender=moderator)
        return html.json_response({'status':'success', 'msg':'user suspended'})
    
    elif action == 'reinstate':
        
        # sanity check, the middleware should disable suspended users loggin in again
        assert moderator != target, 'You may reinstate yourself'

        target.profile.suspended = False
        target.profile.save()
        text = notegen.reinstate(target)
        models.Note.send(target=moderator, content=text, sender=moderator)
        return html.json_response({'status':'success', 'msg':'user reinstated'})
    
    return html.json_response({'status':'error', 'msg':'Invalid action %s' % action})
    
@login_required(redirect_field_name='/openid/login/')
def preview(request):
    "This runs the markdown preview functionality"
    content = request.POST.get('content','no input')[:5000]

    try:
        output = html.generate(content)
    except KeyError, exc:
        # return more userfriendly errors, used for debugging
        output = 'Error: %s' % str(exc)

    return HttpResponse(output, mimetype='text/plain')
