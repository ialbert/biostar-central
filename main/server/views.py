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
        query = query.filter(post_type=POST_QUESTION)    
    if user:
        query = query.filter(author=user)

    query = query.select_related('author','author__profile')

    return query


def index(request):
    "Main page"

    # these are the posts that match the query words
    pids = action.search(request)

    # eventually we will need to order by relevance
    qs = get_posts(request)
    #qs = qs.filter(id__in=pids)
    qs = qs.order_by('-touch_date')
    page  = get_page(request, qs, per_page=20)
    return html.template( request, name='index.html', page=page)

def user_profile(request, uid):
    "User's profile page"
    user = models.User.objects.get(id=uid)
    profile = models.UserProfile.objects.get(user=user)
    profile.writeable = profile.authorize(request.user)
    questions = models.Post.objects.filter(author=user, post_type=POST_QUESTION).select_related('author','author__profile')
    questions = questions.order_by('-score')[:15]
    answers   = models.Post.objects.filter(author=user, post_type=POST_ANSWER).select_related('author', 'author_profile', 'parent__author','parent__author__profile')
    answers   = answers.order_by('-score')[:15]
    notes     = models.Note.objects.filter(target=user).select_related('author', 'author__profile', 'root').order_by('-date')[:15]
    
    answer_count = models.Post.objects.filter(author=user, post_type=POST_ANSWER).count()
    question_count = models.Post.objects.filter(author=user, post_type=POST_QUESTION).count()
    comment_count = models.Post.objects.filter(author=user, post_type=POST_COMMENT).count()

    params = html.Params(question_count=question_count, answer_count=answer_count, comment_count=comment_count)
    return html.template(request, name='user.profile.html',
        user=request.user, profile=profile, selected=user,
        questions=questions,
        answers=answers, notes=notes, params=params)

def user_list(request):
    search  = request.GET.get('search','')[:80] # trim for sanity
    if search:
        query = Q(profile__display_name__icontains=search)
        users = models.User.objects.filter(query).select_related('profile').order_by("-profile__score")
    else:
        users = models.User.objects.select_related('profile').order_by("-profile__score")
    page  = get_page(request, users, per_page=35)
    return html.template(request, name='user.list.html', page=page, rows=7, search=search)

def tag_list(request):
    tags = models.Tag.objects.all().order_by('-count')
    page = get_page(request, tags, per_page=50)
    return html.template(request, name='tag.list.html', page=page)

def badge_list(request):
    badges = models.Badge.objects.filter(secret=False).order_by('name')
    return html.template(request, name='badge.list.html', badges=badges)

def search(request):
    return html.template(request, name='todo.html')

def question_list(request):
    "Lists all the questions" 
    qs = get_questions().filter(answer_count=0)
    page = get_page(request, qs) 
    return html.template(request, name='post.list.html', page=page)

def question_tagged(request, tag_name):
    qs = get_questions().filter(tag_set__name=tag_name)
    page = get_page(request, qs) 
    return html.template(request, name='post.list.html', page=page)

def question_unanswered(request):
    qs = get_questions().filter(answer_count=0)
    page = get_page(request, qs) 
    return html.template(request, name='question.list.html', page=page)
    
def post_show(request, pid):
    "Returns a question with all answers"
    
    qs = models.Post.objects
    question = qs.select_related('children', 'votes').get(id=pid)
        
    #qs = models.Post.all_objects if 'view_deleted' in request.permissions else models.Post.objects
    answers = models.Post.objects.filter(parent=question, post_type=POST_ANSWER).select_related('author', 'author__profile') 
    answers = answers.order_by('-answer_accepted','-score')

    if request.user.is_authenticated():
        notes = models.Note.objects.filter(target=request.user, root=question).all().delete()
        votes = models.Vote.objects.filter(author=request.user, post__id__in=[question.id] + [a.id for a in answers]) 
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
    1/0
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
    assert post_type in POST_REV_MAP, 'Invalid post_type %s' % post_type
  
   
    # this is a readable form of the post type
    post_type_name = POST_REV_MAP.get(post_type,'')

    # select the form factory from the post types
    use_post_form = (post_type_name in POST_FULL_FORM)

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
            text = CREATE_NOTICE.get(post_type, "?")
            post.notify(user=user, text=text)
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
    
    post_type_name = POST_REV_MAP.get(post_type,'')
   
    # select the form factory from the post types
    use_post_form = (post_type_name in POST_FULL_FORM)
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
        text = EDIT_NOTICE.get(post_type, "?")
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

def moderate(request):
    if request.method == 'POST':
        user = request.user
        post_id = int(request.POST.get('post'))
        post = models.Post.objects.get(id=post_id)
        if not post.authorize(user=user, strict=False):
            return html.json_response({'status':'error', 'msg':'You do not have permission to moderate posts.'})        
            
        action = request.POST.get('action')
        action_map = {'close':models.REV_CLOSE, 'reopen':models.REV_REOPEN,
                      'delete':models.REV_DELETE, 'undelete':models.REV_UNDELETE}
        post.moderator_action(action_map[action], user)
        
        return html.json_response({'status':'success', 'msg':'%s performed' % action})
        
    return html.json_response({'status':'error', 'msg':'POST method must be used'})
        


@login_required(redirect_field_name='/openid/login/')
def preview(request):
    "This runs the markdown preview functionality"
    content = request.POST.get('content','no input')

    try:
        output = html.generate(content)
    except KeyError, exc:
        # return more userfriendly errors, used for debugging
        output = 'Error: %s' % str(exc)

    return HttpResponse(output, mimetype='text/plain')
