"""
Biostar views
"""
from main.server import html, models, const, formdef
from main.server.html import get_page

from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login, logout
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q

# the openid association model
from django_openid_auth.models import UserOpenID

# import all constants
from main.server.const import *

def get_questions():
    "Returns a common queryset that can be used to select questions"
    #return models.Question.objects.select_related('post', 'post__author','post__author__profile')
    return models.Post.objects.filter(post_type=POST_QUESTION).select_related('author','author__profile')

def index(request):
    "Main page"
    
    # eventually we will need to order by relevance
    qs = get_questions()
    qs = qs.order_by('-touch_date')
    page  = get_page(request, qs, per_page=20)
    return html.template( request, name='index.html', page=page)

def admin_password_override(request):
    """
    This view is active only if the ALLOW_ADMIN_OVERRIDE setting is True.
    Allows anyone to log in and personify any other user as long as they enter
    the content of SECRET_KEY as password. This needs to be True during testing.
    """
    if request.method=='GET':
        return html.template( request, name='admin.password.override.html')
    else:
        uid = request.POST.get('uid', '')
        passwd = request.POST.get('password', '')
        if passwd == settings.SECRET_KEY:
            user = models.User.objects.get(pk=uid)
            if 'set_type' in request.POST: # Allow setting of user type for testing
                user.profile.type = int(request.POST.get('set_type'))
                user.profile.save()
            user.backend = 'django.contrib.auth.backends.ModelBackend'
            login(request, user)
            return html.redirect('/') 
        
    raise Exception('Invalid login')
        
def user_profile(request, uid):
    "User's profile page"
    user = models.User.objects.get(id=uid)
    profile = models.UserProfile.objects.get(user=user)
    questions = models.Question.objects.filter(post__author=user).select_related('post','post__author','post__author__profile')
    answers = models.Answer.objects.filter(post__author=user).select_related('post','question','question__post','question__post__author','question__post__author__profile')

    return html.template(request, name='user.profile.html',
      user=request.user, profile=profile, selected=user,
      questions=questions.order_by('-post__score'),
      answers=answers.order_by('-post__score'))

def user_list(request):
    search  = request.GET.get('search','')[:80] # trim for sanity
    if search:
        query = Q(first_name__icontains=search) | Q(last_name__icontains=search)
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
    #q.answer_count
    
    qs = get_questions().filter(answer_count=0)
    page = get_page(request, qs) 
    return html.template(request, name='question.list.html', page=page)

def question_tagged(request, tag_name):
    qs = get_questions().filter(post__tag_set__name=tag_name)
    page = get_page(request, qs) 
    return html.template(request, name='question.list.html', page=page)

def question_unanswered(request):
    qs = get_questions().filter(answer_count=0)
    page = get_page(request, qs) 
    return html.template(request, name='question.list.html', page=page)
    
def post_show(request, pid):
    "Returns a question with all answers"
    
    qs = models.Post.objects
    question = qs.filter(post_type=POST_QUESTION).select_related('children').get(id=pid)
    if request.user.is_authenticated():
        question.views += 1
        question.save()
        
    #qs = models.Post.all_objects if 'view_deleted' in request.permissions else models.Post.objects
    qs = models.Post.objects
    answers = qs.filter(parent=question).select_related('author','author__profile','children')
    answers = answers.order_by('-answer_accepted','-score')
    return html.template( request, name='post.show.html', question=question, answers=answers )

def form_revision(post, form):
    "Creates a revision from a form post"
    # sanity check
    assert form.is_valid(), 'form is not valid'
    title   = form.cleaned_data['title']
    content = form.cleaned_data['content']
    tag_string = html.tag_strip(form.cleaned_data['tags'])    
    post.create_revision(content=content, tag_string=tag_string, title=title)
    
@login_required(redirect_field_name='/openid/login/')
def post_parent(request, pid=0):
    """
    Handles parent post related tasks
    """
    
    newpost   = (pid == 0)
    form_data = (request.method == 'POST')
    
    # get post_type 
    post_type = int(request.REQUEST.get('post_type', POST_QUESTION))
    assert post_type in POST_REV_MAP, 'Invalid post_type %s' % post_type
    
    # we have incoming form data for posts with no parents
    if form_data:
        form = formdef.PostForm(request.POST)
        if not form.is_valid():
            return html.template( request, name='edit.post.html', form=form)
        if newpost:
            with transaction.commit_on_success():
                post = models.Post.objects.create(author=request.user, post_type=post_type)
                form_revision(post=post, form=form)
        else:
            post = models.Post.objects.get(pk=pid)
            post.authorize(request)
            form_revision(post=post, form=form)
        return html.redirect('/post/show/%s/%s/' % (post.id, post.slug))
    else:
        if newpost:
            form = formdef.PostForm()
        else:
            post = models.Post.objects.get(pk=pid)
            post.authorize(request)            
            form = formdef.PostForm(initial=dict(title=post.title, content=post.content, tags=post.tag_string))
        return html.template( request, name='edit.post.html', form=form)

@login_required(redirect_field_name='/openid/login/')
def post_child(request, parentid=0, pid=0):
    "Handles actions for posts with parents"

    newpost   = (pid == 0)
    form_data = (request.method == 'POST')
    
    # get post_type 
    post_type = int(request.REQUEST.get('post_type', POST_ANSWER))
    assert post_type in POST_REV_MAP, 'Invalid post_type %s' % post_type
    
    if form_data:
        form = formdef.SimpleForm(request.POST)
        if not form.is_valid():
            return html.template( request, name='edit.post.html', form=form)
            
    
@login_required(redirect_field_name='/openid/login/')
def answer_edit(request, qid, aid=0):
    "Handles answers, requires a question id and answer id"
    
    # get the question that is to be accessed
    post = models.Post.objects.get(pk=qid)

    # decide whether this is new answer 
    newans = (aid == 0)
    
    # it appears to be an edited answer
    edit = not newans

    if edit and request.method == 'GET':
        # editing an existing answer
        answer = models.Post.objects.get(pk=aid)
        # answer edits require an authorization (original author or moderator)
        answer.authorize(request)
        form = AnswerForm( initial=dict(content=answer.content) )
        return html.template( request, name='edit.answer.html', form=form)

    # only POST remains at this point
    assert request.method == 'POST'
    
    # incoming  data
    form = AnswerForm(request.POST)
    if not form.is_valid():
        # invalid form render the errors
        return html.template( request, name='edit.answer.html', form=form)

    # at this point the data is valid
    content = form.cleaned_data['content']

    if newans:
        # new answer for the question
        post = models.Post.objects.create(author=request.user)
        post.create_revision(content=content)
        answer = models.Answer.objects.create(post=post, question=question)
    else:
        # update the answer
        answer = models.Answer.objects.get(pk=aid)
        answer.authorize(request)
        answer.post.create_revision(content=content, author=request.user)

    return html.redirect('/question/show/%s/' % qid)
    
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
def comment_add(request, pid):
    
    parent = models.Post.objects.get(pk=pid)
    content = request.POST['text']
    post = models.Post(author=request.user)
    post.html = post.content = html.sanitize(content, allowed_tags='')
    post.save()
    comment = models.Comment(parent=parent, post=post)
    comment.save()

    try:
        return html.redirect('/question/show/%s/' % (parent.question.id))
    except models.Question.DoesNotExist:
        return html.redirect('/question/show/%s/#%s' % (parent.answer.question.id, parent.answer.id))

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
        author = request.user
        if 'moderate_post' not in request.permissions: # Need to also check for actual mod permissions
            return html.json_response({'status':'error', 'msg':'You do not have permission to moderate posts.'})        

        post_id = int(request.POST.get('post'))
        post = models.Post.objects.get(id=post_id)
        
        action = request.POST.get('action')
        action_map = {'close':models.REV_CLOSE, 'reopen':models.REV_REOPEN,
                      'delete':models.REV_DELETE, 'undelete':models.REV_UNDELETE}
        post.moderator_action(action_map[action], author)
        
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
