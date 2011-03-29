"""
Biostar views
"""
import html
from biostar.server import models
from django import forms
from django.contrib.auth.decorators import login_required
from biostar.libs import postmarkup

def index(request):
    "Main page"
    
    # shows both the 5 freshest and 5 oldest questions 
    # (this is for debugging)
    top = models.Question.objects.all().order_by('-lastedit_date')[:5]
    bot = models.Question.objects.all().order_by('lastedit_date')[:5]
    questions = list(top) + list(bot)
    return html.template( request, name='index.html', questions=questions)

def user(request, uid):
    "User's profile page"
    user = models.User.objects.get(id=uid)
    return html.template(request, name='user.html', selected_user=user)

def users(request):
    users = models.User.objects.all()
    return html.template(request, name='users.html', users=users)

def question_show(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    if request.user.is_authenticated():
        question.post.views += 1
        question.post.save()
    answers  = models.Answer.objects.filter(question=question).select_related()

    return html.template( request, name='question.html', question=question, answers=answers )

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
        form = QuestionForm(initial=dict(title=question.title, content=question.post.bbcode, tags=tags))
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
        post = models.Post.objects.create(bbcode=content, author=request.user)
        question = models.Question.objects.create(post=post, title=title)
        question.tags.add(*tags)
    else:
        # editing existing question
        question = models.Question.objects.get(pk=pid)
        post = question.post
        question.title, post.bbcode, post.lastedit_user = title, content, request.user
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
    
    if aid != 0:
        # need to authenticate user access to this answer
        # no authentication for now
        answer = models.Answer.objects.get(pk=aid)
       
    if request.method == 'GET':
        # editing an existing answer
        form = AnswerForm( initial=dict(content=answer.post.bbcode) )
        return html.template( request, name='edit.answer.html', form=form)

    elif request.method == 'POST':
        # incoming form data
        form = AnswerForm(request.POST)

        if form.is_valid():
            content = form.cleaned_data['content']
            if aid==0:
                # new answer for the question
                post = models.Post(bbcode=content, author=request.user)
                post.save()
                answer = models.Answer(post=post, question=question)
                answer.save()
            else:
                # the answer is already generated in this view
                answer.post.bbcode = content
                answer.post.save()
            # go to the question
            return html.redirect('/question/show/%s/' % qid)
        else:
            # invalid form submission, render the errors
            return html.template( request, name='edit.answer.html', form=form)

def vote(request):
    "Handles all voting on posts"
    if request.method == 'POST':
        
        author = request.user
        if not author.is_authenticated():
            return html.json_response({'status':'error', 'msg':'You must be logged in to vote'})
        
        post_id = int(request.POST.get('post'))
        post = models.Post.objects.get(id=post_id)
        
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
        
