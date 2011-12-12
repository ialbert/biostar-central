"""
Biostar views
"""
import difflib
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

def index(request):
    "Main page"
    
    params = html.Params()
    params.parse(request)
    
    logger.info('testing the logger')

    if params.q:
        pids  = action.search(params.q)
        posts = get_posts(request).filter(id__in=pids).order_by('-touch_date')
    else:        
        # no search was performed, get the latest questions
        posts = get_posts(request).filter(type=POST_QUESTION).order_by('-touch_date')
        
    page  = get_page(request, posts, per_page=20)
    return html.template(request, name='index.html', page=page, params=params)

def post_list_filter(request, uid=0, word=None):
    post_type = {  'questions': POST_QUESTION, 'answers':POST_ANSWER, 'comments': POST_COMMENT }.get(word)
    return post_list(request, uid=uid, post_type=post_type)

def post_list(request, uid=0, post_type=None):
    params = html.Params()

    posts = get_posts(request).filter(type=post_type)
    if uid:
        user = models.User.objects.filter(id=uid).select_related('profile').all()[0]
        posts = posts.filter(author=user)
        params.setr('Filter: %s' % user.profile.display_name)
    posts = posts.order_by('-lastedit_date')
    page  = get_page(request, posts, per_page=20)
    return html.template( request, name='post.list.html', page=page, params=params)

def user_profile(request, uid):
    "User's profile page"

    user = request.user
    target  = models.User.objects.get(id=uid)
    
    if target == user:
        notes = models.Note.objects.filter(target=target).select_related('author', 'author__profile', 'root').order_by('-date')
        page  = get_page(request, notes, per_page=20)
        # we evalute it here so that subsequent status updates won't interfere
        page.object_list = list(page.object_list)
        models.Note.objects.filter(target=target).update(unread=False)
        models.UserProfile.objects.filter(user=target).update(new_messages=0)
    else:
        page = None

    # we need to collate and count the awards
    awards = models.Award.objects.filter(user=target).select_related('badge').order_by('-date')
    
    answer_count = models.Post.objects.filter(author=target, type=POST_ANSWER).count()
    question_count = models.Post.objects.filter(author=target, type=POST_QUESTION).count()
    comment_count = models.Post.objects.filter(author=target, type=POST_COMMENT).count()
    post_count = models.Post.objects.filter(author=target).count()
    vote_count = models.Vote.objects.filter(author=target).count()
    award_count = models.Award.objects.filter(user=target).count()
    
    params = html.Params(question_count=question_count, answer_count=answer_count, 
        comment_count=comment_count, post_count=post_count, vote_count=vote_count, award_count=award_count)
    
    # some information is only visible to the user
    target.writeable = auth.authorize_user_edit(target=target, user=user, strict=False)
    target.showall  = (target == user)
    
    return html.template(request, name='user.profile.html', awards=awards,
        user=request.user,target=target, params=params, page=page)

def user_list(request):
    search  = request.GET.get('m','')[:80] # trim for sanity
    if search:
        query = Q(profile__display_name__icontains=search)
        users = models.User.objects.filter(query).select_related('profile').order_by("-profile__score")
    else:
        users = models.User.objects.select_related('profile').order_by("-profile__score")
    page  = get_page(request, users, per_page=24)
    return html.template(request, name='user.list.html', page=page)

def tag_list(request):
    tags = models.Tag.objects.all().order_by('-count')
    page = get_page(request, tags, per_page=50)
    return html.template(request, name='tag.list.html', page=page)

def badge_list(request):
    badges = models.Badge.objects.filter(secret=False).order_by('-count', '-type')
    return html.template(request, name='badge.list.html', badges=badges)

def question_unanswered(request, uid=0, post_type=None):
    "Lists all the questions"
    params = html.Params()
    params.setr('Filter: unanswered')
    qs = get_posts(request).filter(answer_count=0)
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
    
    qs = models.Post.objects
    question = qs.select_related('children', 'votes').get(id=pid)
            
    #qs = models.Post.all_objects if 'view_deleted' in request.permissions else models.Post.objects
    answers = models.Post.objects.filter(parent=question, type=POST_ANSWER).select_related('author', 'author__profile') 
    answers = list(answers.order_by('-accepted','-score'))
    
    # add the writeable attribute to each post
    all = [ question ] + answers
    for post in all:
        auth.authorize_post_edit(post=post, user=request.user, strict=False)
    
    if request.user.is_authenticated():
        notes = models.Note.objects.filter(target=request.user, post=question).all().delete()
        votes = models.Vote.objects.filter(author=request.user, post__id__in=[ question.id ] + [a.id for a in answers] ) 
        
        # updates the viewcounter once within a session, Alex says to move to IP based counting TODO
        viewed = request.session.get(VIEWED_KEY, set())
        if question.id not in viewed:
            viewed.add(question.id)
            question.views += 1
            question.save()
            request.session[VIEWED_KEY] = viewed
    else:
        notes, votes = [], []
        
        
    up_votes = set(vote.post.id for vote in votes if vote.type == const.VOTE_UP)
    down_votes = set(vote.post.id for vote in votes if vote.type == const.VOTE_DOWN)
    
    #return html.template( request, name='post.html', question=question, answers=answers, up_votes=up_votes, down_votes=down_votes )
 
    return html.template( request, name='post.show.html', question=question, answers=answers, up_votes=up_votes, down_votes=down_votes )

def post_redirect(post, anchor=None):
    """
    Shows a post in full context
    """
    # get the root of a post
    pid, slug = post.root.id, post.root.slug
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
    use_post_form = (post_type not in POST_CONTENT_ONLY)

    if use_post_form:
        factory = formdef.PostForm
    else:
        factory = formdef.ContentForm
    
    # find the parent if it exists
    if parentid:
        parent = models.Post.objects.get(pk=parentid)
    else:
        parent = None

    # deal with new post creation first
    if newpost:
        # this here is to customize the output
        params = html.Params(title="New post", use_post_form=use_post_form ) 

        if form_data:
            form = factory(request.POST)
            if not form.is_valid():
                return html.template(request, name='post.edit.html', form=form, params=params)
            params = dict(author=user, type=post_type, parent=parent, root=parent, creation_date=datetime.now())
            params.update(form.cleaned_data)            
            with transaction.commit_on_success():
                post = models.Post.objects.create(**params)
                post.set_tags()
                models.create_revision(post)
            return post_redirect(post)
        else:
            form = factory()
            return html.template( request, name='post.edit.html', form=form, params=params)

    #
    # at this point we are dealing with a post editing action
    #
    assert pid, 'Only post modification should follow after this point'
    
    # when we edit a post we keep the original post type (for now)
    post = models.Post.objects.get(pk=pid)
    parent = post.parent
    post_type = post.type
   
    # select the form factory from the post types
    use_post_form = (post_type not in POST_CONTENT_ONLY)
    if use_post_form:
        factory = formdef.PostForm
    else:
        factory = formdef.ContentForm

    # verify that this user may indeed modify the post
    auth.authorize_post_edit(post=post, user=request.user, strict=True)
    
    params = html.Params(title="Edit %s" % post.get_type_display(), use_post_form=use_post_form) 

    # no form data coming, return the editing form
    if not form_data:
        form = factory(initial=dict(title=post.title, content=post.content, tag_val=post.tag_val))
        return html.template(request, name='post.edit.html', form=form, params=params)
    else:
        # we have incoming form data for posts
        form = factory(request.POST)
        if not form.is_valid():
            return html.template(request, name='post.edit.html', form=form, params=params)
         
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
    
def vote(request):
    "Handles all voting on posts"
    if request.method == 'POST':
        
        author = request.user
        if not author.is_authenticated():
            return html.json_response({'status':'error', 'msg':'You must be logged in to vote'})
        
        post_id = int(request.POST.get('post'))
        post = models.Post.objects.get(id=post_id)
        
        if post.author == author:
            return html.json_response({'status':'error', 'msg':'You cannot vote on your own post'})
        
        type = int(request.POST.get('type'))
        
        old_vote = post.get_vote(author, type)
        
        if old_vote:
            old_vote.delete()
            return html.json_response({
                'status':'success',
                'msg':'%s removed' % old_vote.get_type_display()})
        else:
            vote = post.add_vote(author, type)
            if type in models.OPPOSING_VOTES: 
                # Remove an opposing vote if it exists
                post.remove_vote(author, models.OPPOSING_VOTES[type])
            return html.json_response({
                'status':'success',
                'msg':'%s added' % vote.get_type_display()})
                    

    return html.json_response({'status':'error', 'msg':'POST method must be used'})

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
