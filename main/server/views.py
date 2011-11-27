"""
Biostar views
"""
from main.server import html, models, const, formdef, action
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


def get_posts(request, post_type=POST_QUESTION, user=None):
    "Returns a common queryset that can be used to select questions"
    
    query = models.Post.objects
    
    if post_type:
        query = query.filter(post_type=post_type)  

    if user:
        query = query.filter(author=user)

    query = query.select_related('author','author__profile')

    return query

def index(request):
    "Main page"

    # this will contain the query if it was sent 
    query = request.REQUEST.get('q','')
    pids  = action.search(query)
    params = html.Params()
    if query:
        params.remind = 'Searching for: %s' % query
        qs = get_posts(request, post_type=None)
        qs = qs.filter(id__in=pids)
    else:
        qs = get_posts(request)
        
    qs  = qs.order_by('-touch_date')
    page  = get_page(request, qs, per_page=20)
    return html.template( request, name='index.html', page=page, params=params)

def post_list_filter(request, uid=0, word=None):
    post_type = {  'questions': POST_QUESTION, 'answers':POST_ANSWER, 'comments': POST_COMMENT }.get(word)
    return post_list(request, uid=uid, post_type=post_type)

def post_list(request, uid=0, post_type=None):
    params = html.Params()

    posts = get_posts(request, post_type=post_type)
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
        page  = get_page(request, notes, per_page=10)
        # we evalute it here so that subsequent status updates won't interfere
        page.object_list = list(page.object_list)
        models.Note.objects.filter(target=target).update(unread=False)
    else:
        page = None

    # we need to collate and count the awards
    awards = models.Award.objects.filter(user=target).select_related('badge').order_by('-date')
    
    answer_count = models.Post.objects.filter(author=target, post_type=POST_ANSWER).count()
    question_count = models.Post.objects.filter(author=target, post_type=POST_QUESTION).count()
    comment_count = models.Post.objects.filter(author=target, post_type=POST_COMMENT).count()
    post_count = models.Post.objects.filter(author=target).count()
    vote_count = models.Vote.objects.filter(author=target).count()
    award_count = models.Award.objects.filter(user=target).count()
    
    params = html.Params(question_count=question_count, answer_count=answer_count, 
        comment_count=comment_count, post_count=post_count, vote_count=vote_count, award_count=award_count)
    
    # it the target personal information writable
    target.editable   = target.profile.editable(user)
    target.authorized = target.profile.authorize(user)
    return html.template(request, name='user.profile.html',
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

def search(request):
    return html.template(request, name='todo.html')

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
    answers = models.Post.objects.filter(parent=question, post_type=POST_ANSWER).select_related('author', 'author__profile') 
    answers = answers.order_by('-answer_accepted','-score')

    if request.user.is_authenticated():
        notes = models.Note.objects.filter(target=request.user, post=question).all().delete()
        votes = models.Vote.objects.filter(author=request.user, post__id__in=[ question.id ] + [a.id for a in answers] ) 
        question.views += 1
        question.save()
    else:
        notes, votes = [], []
        
        
    up_votes = set(vote.post.id for vote in votes if vote.type == const.VOTE_UP)
    down_votes = set(vote.post.id for vote in votes if vote.type == const.VOTE_DOWN)
     
    return html.template( request, name='post.show.html', question=question, answers=answers, up_votes=up_votes, down_votes=down_votes )

def show_post(post, anchor=None):
    """
    Shows a post in full context
    """
    # get the root of a post
    root = post.get_root()
    pid, slug = root.id, root.slug
    anchor = anchor or post.id
    url = '/post/show/%s/%s/#%s' % (pid, slug, anchor)
    return html.redirect(url)

def process_form(post, form, user):
    "Creates a revision from a form post"
    # sanity check
    assert form.is_valid(), 'form is not valid'
    title   = form.cleaned_data.get('title','')
    content = form.cleaned_data.get('content', '')
    tag_string = form.cleaned_data.get('tags_string', '')
    tag_string = html.tag_strip(tag_string)   
    post.create_revision(content=content, tag_string=tag_string, title=title, author=user)
    
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
  
    # this is a readable form of the post type
    post_type_name = POST_MAP.get(post_type,'')

    # select the form factory from the post types
    use_post_form = (post_type in POST_FULL_FORM)

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
        params = html.Params(title="New %s" % post_type_name, use_post_form=use_post_form ) 

        if form_data:
            form = factory(request.POST)
            if not form.is_valid():
                return html.template(request, name='post.edit.html', form=form, params=params)
            params = dict(author=user, post_type=post_type, parent=parent, creation_date=datetime.now())
            params.update(form.cleaned_data)            
            post = models.Post.objects.create(**params)
            process_form(post, form, user=request.user)
            return show_post(post)

        else:
            form = factory()
            return html.template( request, name='post.edit.html', form=form, params=params)

    # at this point we are dealing with a post editing aciont
    assert pid, 'Only post modification should follow after this point'
    
    # when we edit a post we keep the original post type (for now)
    post = models.Post.objects.get(pk=pid)
    parent = post.parent
    post_type = post.post_type
    
    post_type_name = POST_MAP.get(post_type,'')
   
    # select the form factory from the post types
    use_post_form = (post_type in POST_FULL_FORM)
    if use_post_form:
        factory = formdef.PostForm
    else:
        factory = formdef.ContentForm

    # verify that this user may indeed modify the post
    post.authorize(user=request.user, strict=True)
    
    params = html.Params(title="Edit %s" % post_type_name, use_post_form=use_post_form) 

    # no form data coming, return the editing form
    if not form_data:
        form = factory(initial=dict(title=post.title, content=post.content, tag_string=post.tag_string))
        return html.template(request, name='post.edit.html', form=form, params=params)
    else:
        # we have incoming form data for posts
        form = factory(request.POST)
        if not form.is_valid():
            return html.template(request, name='post.edit.html', form=form, params=params)
        process_form(post=post, form=form, user=request.user)        
        return show_post(post)    

def revision_list(request, pid):
    post = models.Post.objects.get(pk=pid)
    revisions = list(post.revisions.order_by('date')) # Oldest first, will need to be reversed later
    
    # We need to annotate the revisions with exactly what was changed
    # 'verb' is a word for the action box to describe the revision
    # 'modified' is a list (title, content, tag_string) of boolean values for if it was changed
    def revision_data(rev):
        return rev.title, rev.content, rev.tag_string
    last_data = revision_data(revisions[0])
    revisions[0].verb = 'created'
    revisions[0].modified = [True, True, True] # Always display the first revision
    for revision in revisions[1:]:
        if revision.action:
            revision.verb = 'performed'
        else:
            revision.verb = 'edited'
            data = revision_data(revision)
            revision.modified = [current != last for current, last in zip(data, last_data)]
            last_data = data
    revisions.reverse()
    
    return html.template(request, name='revision.list.html', revisions=revisions, post=post)
   
@login_required(redirect_field_name='/openid/login/')
def add_comment(request, pid):
    "Adds a comment"
    parent  = models.Post.objects.get(pk=pid)
    content = request.POST['text'].strip()
    if len(content)<1:
        messages.warning(request, 'Comment too short!')
        return show_post(parent)
        
    comment = models.Post(author=request.user, parent=parent, post_type=POST_COMMENT, creation_date=datetime.now())
    comment.save()
    comment.create_revision(content=content)
    
    return show_post(comment)
    
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
    
    if request.method != 'POST':
        return html.json_response({'status':'error', 'msg':'Only POST requests are allowed'})        
   
    moderator = request.user
    post = models.Post.objects.get(id=pid)
    if not post.authorize(user=moderator, strict=False):
        return html.json_response({'status':'error', 'msg':'You do not have permission to moderate this post.'})        
        
    action_map = {'close':models.REV_CLOSE, 'reopen':models.REV_REOPEN,
                  'delete':models.REV_DELETE, 'undelete':models.REV_UNDELETE}
    
    post.moderator_action(action_map[action], moderator)
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
        text = 'suspended %s' % target.profile.display_name
        models.Note.send(target=target, content=text, sender=moderator)
        return html.json_response({'status':'success', 'msg':'user suspended'})
    
    elif action == 'reinstate':
        
        # sanity check, the middleware should disable suspended users loggin in again
        assert moderator != target, 'You may reinstate yourself'

        target.profile.suspended = False
        target.profile.save()
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
