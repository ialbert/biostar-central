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

    return html.template( request, name='index.html', params=params )

def question(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    answers  = models.Answer.objects.filter(question=question).select_related()

    params = html.Params(question=question, answers=answers )
    return html.template( request, name='question.html', params=params )

class PostForm(forms.Form):
    content = forms.CharField(max_length=1000)
    parent  = forms.IntegerField(required=False, initial=0)

def newpost(request):
    "Handles all new posts: questions, answers and comments"
    if request.method == 'POST':
        form = PostForm(request.POST)
        
        if form.is_valid(): # All validation rules pass
            parent  = form.cleaned_data['parent']
            content = form.cleaned_data['content']

            # create the HTML from the bbcode
            parse = postmarkup.create()
            body  = parse(content)
         
            # in this demo all new posts go under user number 1
            author = models.User.objects.get(id=1)
            
            post = models.Post(bbcode=content, html=body, author=author, lastedit_user=author)
            post.save()
            
            if not parent:
                # no parent means new question
                question = models.Answer(post=post)
                question.save()
            else:
                # an answer to a question
                question = models.Question.objects.get(id=parent)
                answer = models.Answer(question=question, post=post)
                answer.save()
            return html.redirect('/question/%s/' % question.id) 
    else:
        return html.redirect('/about/') 
