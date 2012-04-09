"""
Generating awards
"""
    
from django.conf import settings
from main.server import models, html, notegen
from main.server.const import *
from django.contrib import messages

#models.Award.objects.all().delete()

def init_badges():
    "Initializes badges "
    models.Badge.get_or_create(name='Teacher', description='')
    pass

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
        
    pairs = [
        ('Teacher', models.Post.objects.filter(author=user, score__gt=0).count),
        ('Supporter', models.Vote.objects.filter(author=user).count),
        ('Nice Question', models.Post.objects.filter(author=user, score__gt=10).count),
        ('Famous Question', models.Post.objects.filter(author=user, score__gt=250).count),
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