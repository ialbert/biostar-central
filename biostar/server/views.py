"""
Biostar views
"""
import html
from biostar.server import models
from django import forms
from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login
from django.conf import settings
from taggit.models import Tag

from django.http import HttpResponse
import markdown
from django.views.decorators.csrf import csrf_exempt


def index(request):
    "Main page"
    
    if request.user.is_authenticated():
        merge_accounts(request)

    # shows both the 5 freshest and 5 oldest questions 
    # (this is for debugging)
    top = models.Question.objects.all().order_by('-lastedit_date')[:5]
    bot = models.Question.objects.all().order_by('lastedit_date')[:5]
    questions = list(top) + list(bot)
    return html.template( request, name='index.html', questions=questions)

def test_login(request):
    "Allows for an automatic login used during testing"
    if settings.DEBUG:
        user = authenticate(username='testuser', password='test$123')
        login(request, user)
    else:
        raise Exception('Invalid login')

@transaction.commit_on_success
def merge_accounts(request):
    "Attempts to merge user accounts if emails match"
    users = list(models.User.objects.filter(email=request.user.email))
    if len(users)>1:
        source, target = users[0], users[-1]
        
        models.Post.objects.filter(author=source).update(author=target)
        models.Vote.objects.filter(author=source).update(author=target)
        
        # needs a one step transfer of all attributes
        p1 = models.UserProfile.objects.get(user=source)
        p2 = models.UserProfile.objects.get(user=target)
        p2.score = p1.score
        p2.save()
        
        # disable the old user
        source.set_unusable_password()

def user_profile(request, uid):
    "User's profile page"
    user = models.User.objects.get(id=uid)
    questions = models.Question.objects.filter(post__author=user)

    return html.template(request, name='user.profile.html', selected_user=user, questions=questions)

def get_page(request, obj_list, per_page=25):
    "A generic paginator"

    paginator = Paginator(obj_list, per_page) 
    try:
        pid = int(request.GET.get('page', '1'))
    except ValueError:
        pid = 1

    try:
        page = paginator.page(pid)
    except (EmptyPage, InvalidPage):
        page = paginator.page(paginator.num_pages)
    
    return page

def user_list(request):
    users = models.User.objects.all()
    page  = get_page(request, users)
    return html.template(request, name='user.list.html', page=page)

def tag_list(request):
    tags = Tag.objects.all()
    page = get_page(request, tags)
    return html.template(request, name='tag.list.html', page=page)

def badge_list(request):
    return html.template(request, name='todo.html')

def question_unanswered(request):
    return html.template(request, name='todo.html')

def search(request):
    return html.template(request, name='todo.html')

def question_list(request):
    "Lists all the questions"
    all  = models.Question.objects.all()
    page = get_page(request, all) 
    return html.template(request, name='question.list.html', page=page)

def question_show(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    if request.user.is_authenticated():
        question.post.views += 1
        question.post.save()
    answers  = models.Answer.objects.filter(question=question).select_related()

    return html.template( request, name='question.show.html', question=question, answers=answers )

# question form and its default values
_Q_TITLE, _Q_CONTENT, _Q_TAG = 'Question title', 'question content', 'tag1'
class QuestionForm(forms.Form):
    "A form representing a new question"
    title   = forms.CharField(max_length=250,  initial=_Q_TITLE)
    content = forms.CharField(max_length=5000, initial=_Q_CONTENT)
    tags    = forms.CharField(max_length=250,  initial=_Q_TAG)
    
    def clean(self):
        "Custom validator for the question"
        if self.cleaned_data['tags'] == _Q_TAG:
            raise forms.ValidationError("Please change the tag from default value")
    
        if self.cleaned_data['content'] == _Q_CONTENT:
            raise forms.ValidationError("Please change content from default value")
        
        if self.cleaned_data['title'] == _Q_TITLE:
            raise forms.ValidationError("Please change title from default value")

        return self.cleaned_data

# editing questions/answers and comments can be merged into a single post handler
# for now these are separete to allow us to identify what each step needs

@login_required(redirect_field_name='/openid/login/')
def question_edit(request, pid=0):
    "Handles questions"
    
    # pid == 0 means a new question
    asknew = pid == 0
    edit = not asknew

    # rather than nesting ifs lets deal with each case separately 
    # then exit as soon as we got enough information
    if asknew:
        params = html.Params(subheader='Ask a question', button='Ask your question')
    else:
        params = html.Params(subheader='Edit question', button='Submit changes')

    # looks like a new question
    if asknew and request.method == 'GET':
        form = QuestionForm()
        return html.template( request, name='edit.question.html', form=form, params=params)
    
    # looks like starting an edit request for an existing question
    if edit and request.method == 'GET':
        question = models.Question.objects.get(pk=pid)
        question.authorize(request)
        tags = " ".join([ tag.name for tag in question.tags.all() ])
        form = QuestionForm(initial=dict(title=question.title, content=question.post.content, tags=tags))
        return html.template( request, name='edit.question.html', form=form, params=params)
  
    # we can only deal with POST after this point
    assert  request.method == 'POST'
    
    # check the form for validity
    form = QuestionForm(request.POST)
    if not form.is_valid():
        return html.template( request, name='edit.question.html', form=form, params=params)

    # at this point the form is valid
    title   = form.cleaned_data['title']
    content = form.cleaned_data['content']
    tags    = form.cleaned_data['tags'].split()
            
    if asknew:
        # generate the new question
        post = models.Post.objects.create(author=request.user)
        post.set(content)
        question = models.Question.objects.create(post=post, title=title)
        question.tags.add(*tags)
    else:
        # editing existing question
        question = models.Question.objects.get(pk=pid)
        post = question.post
        post.set(content)
        question.title, post.lastedit_user = title, request.user
        question.tags.add(*tags)
        question.save(), post.save()

    # show the question
    return html.redirect('/question/show/%s/' % question.id) 
            

# answer/comment form and its default values
class AnswerForm(forms.Form):
    "A form representing an answer or comment"
    content = forms.CharField(max_length=5000)

@login_required(redirect_field_name='/openid/login/')
def answer_edit(request, qid, aid=0):
    "Handles answers, question id and answer id"
    
    # get the question that is to be accessed
    question = models.Question.objects.get(pk=qid)

    newans = aid == 0
    edit = not newans

    if edit and request.method == 'GET':
        # editing an existing answer
        answer = models.Answer.objects.get(pk=aid)
        answer.authorize(request)
        form = AnswerForm( initial=dict(content=answer.post.content) )
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
        post.set(content)
        answer = models.Answer.objects.create(post=post, question=question)
    else:
        # update the answer
        answer = models.Answer.objects.get(pk=aid)
        answer.authorize(request)
        answer.post.set(content)

    return html.redirect('/question/show/%s/' % qid)
   
@login_required(redirect_field_name='/openid/login/')
def comment_add(request, pid):
    
    parent = models.Post.objects.get(pk=pid)
    content = request.POST['text']
    post = models.Post(author=request.user)
    post.html = post.content = html.sanitize(content, allowed_tags='')
    post.save()
    comment = models.Comment(parent=parent, post=post)
    comment.save()

    if parent.question_set.count(): # Post is a question
        return html.redirect('/question/show/%s/' % (parent.question_set.all()[0].id))
    else:
        return html.redirect('/question/show/%s/' % (parent.answer_set.all()[0].question.id))
    
    
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

@csrf_exempt
def markdown_preview(request):
    source_text = request.REQUEST['source_text'] # May need to be sanitized here
    html = markdown.markdown(source_text, safe_mode='remove')
    return HttpResponse(html, mimetype='text/plain')

