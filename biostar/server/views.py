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

_Q_TITLE, _Q_CONTENT, _Q_TAG = 'Question title', 'question content', 'tag1'
class QuestionForm(forms.Form):
    "A form representing a new question"
    title   = forms.CharField(max_length=250,  initial=_Q_TITLE)
    content = forms.CharField(max_length=5000, initial=_Q_CONTENT)
    tags    = forms.CharField(max_length=250,  initial=_Q_TAG)
    
    def clean(self):
        "Custom validator for the form"
        if self.cleaned_data['tags'] == _Q_TAG:
            raise forms.ValidationError("Please change the tag from default value")
    
        if self.cleaned_data['content'] == _Q_CONTENT:
            raise forms.ValidationError("Please change content from default value")
        
        if self.cleaned_data['title'] == _Q_TITLE:
            raise forms.ValidationError("Please change title from default value")

        return self.cleaned_data

class AnswerForm(forms.Form):
    "A form representing an answer or comment"
    parent  = forms.IntegerField(0)
    content = forms.CharField(max_length=5000)

# editing questions/answers and comments can be merged into a single post handler
# for now these are separete to allow us to identify what each step needs

def question_show(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    if request.user.is_authenticated():
        question.post.views += 1
        question.post.save()
    answers  = models.Answer.objects.filter(question=question).select_related()

    return html.template( request, name='question.html', question=question, answers=answers )

@login_required(redirect_field_name='/openid/login/')
def question_edit(request, pid=0):
    "Handles new questions"
    
    if pid==0:
        # looks like a new question
        params = html.Params(subheader='Ask a question', button='Ask your question')
        form = QuestionForm()
    else:
        # looks like editing an existing question
        params = html.Params(subheader='Edit question', button='Submit changes')
        
        # we will need to validate edit access to the question by the author
        # no authorization check for now
        question = models.Question.objects.get(pk=pid)
        tags = " ".join([ tag.name for tag in question.tags.all() ])
        form = QuestionForm(initial=dict(title=question.title, content=question.post.bbcode, tags=tags))
    
    if request.method == 'GET':
        return html.template( request, name='edit.question.html', form=form, params=params)

    elif request.method == 'POST':
        # incoming data posted
        form = QuestionForm(request.POST)
    
        if form.is_valid():
            # generate the new question
            title   = form.cleaned_data['title']
            content = form.cleaned_data['content']
            tags    = form.cleaned_data['tags'].split()
            
            if pid == 0:
                # new question
                post = models.Post(bbcode=content, author=request.user)
                post.save()
                question = models.Question(post=post, title=title)
                question.save()
                question.tags.add(*tags)
            else:
                # editing existing question

                # first it needs to validate access to the question then update
                question = models.Question.objects.get(pk=pid)
                # no authorization check for now
                post = question.post
                question.title, post.bbcode, post.lastedit_user = title, content, request.user
                question.tags.add(*tags)
                question.save(), post.save()

            # redirect to the question
            return html.redirect('/question/%s/show/' % question.id) 
        else:
            # return form with error message
            return html.template( request, name='edit.question.html', form=form, params=params)
    
'''
@login_required(redirect_field_name='/openid/login/')
def newpost(request):
    "Handles all new posts: questions, answers and comments"
    if request.method == 'POST':
        form = PostForm(request.POST)
    else:
        raise HttpError("/")
    if form.is_valid(): # All validation rules pass
        parent  = form.cleaned_data['parent']
        title   = form.cleaned_data['title']
        content = form.cleaned_data['content']
    
        post = models.Post(bbcode=content, author=request.user)
        post.save()
            
            if not parent:
                # no parent means it is a new question
                question = models.Question(post=post, title=title)
                question.save()
            else:
                # an answer to a question
                question = models.Question.objects.get(id=parent)
                answer = models.Answer(question=question, post=post)
                answer.save()
            return html.redirect('/question/%s/' % question.id) 
    else:
        return html.redirect('/about/') 

'''

def vote(request):
    "Handles all voting on posts"
    if request.method == 'POST':
        
                
        # in this demo all votes go under user number 1
        #author = models.User.objects.get(id=1)
        
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
            if type in models.OPPOSING_VOTES: # Remove an opposing vote if it exists
                post.remove_vote(author, models.OPPOSING_VOTES[type])
            return html.json_response({
                'status':'success',
                'msg':'%s added' % vote.get_type_display()})
                    

    return html.json_response({'status':'error', 'msg':'POST method must be used'})
        
