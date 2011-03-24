"""
Biostar views
"""
import html
from biostar.server import models

def index(request):
    "Main page"
    questions = models.Question.objects.all()[:10]
    
    params = html.Params(questions=questions )

    return html.template( request, name='index.html', params=params )

def question(request, pid):
    "Returns a question with all answers"
    question = models.Question.objects.get(id=pid)
    answers  = models.Answer.objects.filter(question=question).select_related()
    
    print answers

    params = html.Params(question=question, answers=answers )

    return html.template( request, name='question.html', params=params )