"""
Biostar views
"""
import html
from biostar.server import models
from django import forms
from biostar.libs import postmarkup

def index(request):
    "Main page"
    questions = models.Question.objects.all()[:10]
    
    params = html.Params(questions=questions )

    return html.template( request, name='index.html', questions=questions)

def question(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    answers  = models.Answer.objects.filter(question=question).select_related()

    params = html.Params(question=question, answers=answers )
    return html.template( request, name='question.html', params=params )

class PostForm(forms.Form):
    content = forms.CharField(max_length=1000)
    title   = forms.CharField(max_length=250, required=False)
    parent  = forms.IntegerField(required=False, initial=0)

def newpost(request):
    "Handles all new posts: questions, answers and comments"
    if request.method == 'POST':
        form = PostForm(request.POST)
        
        if form.is_valid(): # All validation rules pass
            parent  = form.cleaned_data['parent']
            title   = form.cleaned_data['title']
            content = form.cleaned_data['content']

            # in this demo all new posts go under user number 1
            author = models.User.objects.get(id=1)
            
            post = models.Post(bbcode=content, author=author)
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
            return html.json_response({'status':'success', 'msg':'%s removed' % old_vote.get_type_display()})
        else:
            vote = models.Vote(post=post, author=author, type=type)
            vote.save()
            return html.json_response({'status':'success', 'msg':'%s added' % vote.get_type_display()})
                    

    return html.json_response({'status':'error', 'msg':'POST method must be used'})
        
