"""
Generating awards
"""
from datetime import date, timedelta

from django.conf import settings
from main.server import models, html, notegen
from main.server.const import *
from django.contrib import messages

#models.Award.objects.all().delete()

def create(request, user, badge):
    award = models.Award.objects.create(user=user, badge=badge)
    text = notegen.awardnote(award.badge)
    note = models.Note.objects.create(sender=user, target=user, content=text, url=award.badge.get_absolute_url() )
    messages.info(request, note.html)

def instant(request):
    """
    Produces an instant award if applicable and returns.
    These awards may be granted during active sessions
    """
    user = request.user

    badges = dict( [ (b.name, b) for b in  models.Badge.objects.all() ] )
    awards = set( models.Award.objects.filter(user=user).values_list('badge__name', flat=True).distinct() )
    
    def apply_award(name, func):
        "Applies an award"
        badge = badges.get(name)
        if badge and badge.name not in awards and func():
            create(request, user=user, badge=badge)
            return True
        return False
    
    def civic_duty():
        "Selector for Civid Duty badge"
        return models.Vote.objects.filter(author=user, type=VOTE_UP).count() > 300
    
    def commentator():
        return models.Post.objects.filter(author=user,type=POST_COMMENT).count() > 10
    
    def guru():
        return models.Post.objects.filter(author=user, type=POST_ANSWER, score__gt=40).count()
        
    def student():
        return models.Post.objects.filter(author=user, type=POST_QUESTION, score__gt=0).count()
    
    def teacher():
        return models.Post.objects.filter(author=user, type=POST_ANSWER, score__gt=0).count()
        
    def editor():
        return models.PostRevision.objects.filter(author=user).count()
        
    def supporter():
        return models.Vote.objects.filter(author=user).count()
        
    def nice_question():
        return models.Post.objects.filter(author=user, score__gt=10).count()
    
    def pundit():
        return models.Post.objects.filter(author=user, type=POST_COMMENT, score__gt=10).count()
      
    def yearling():
        start = date.today() - timedelta(days=365)
        return models.User.objects.filter(id=user.id, date_joined__lt=start, profile__score__gt=20).count()
    
    def autobiographer():
        p = user.profile
        
        return p.about_me and p.website and p.location
       
       
    pairs = [
        ('Teacher', teacher),
        ('Guru', guru),
        ('Autobiographer', autobiographer),
        ('Yearling', yearling),
        ('Pundit', pundit),
        ('Editor', editor),
        ('Student', student),
        ('Commentator', commentator),
        ('Supporter', supporter),
        ('Nice Question', nice_question ),
        ('Civic Duty', civic_duty),
        ]
    
    for name, func in pairs:
        if apply_award(name, func):
            return
        
    # Famous question has its own verification as it is a multi award
    badge = badges.get('Famous Question')
    if badge:
        badge_count = models.Award.objects.filter(user=user, badge=badge).count()
        post_count  = models.Post.objects.filter(author=user, views__gt=5000).count()
        if badge_count < post_count:
            create(request, user=user, badge=badge)
            return
        