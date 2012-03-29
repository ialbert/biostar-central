"""
Generating awards
"""
    
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

def check(badge, awards, func):
    return badge and badge.name not in awards and func()
        
def instant(request):
    "Produces an instant award if applicable and returns"
    user = request.user

    badges = dict( [ (b.name, b) for b in  models.Badge.objects.all() ] )
    awards = set( models.Award.objects.values_list('badge__name', flat=True).distinct() )
    
    badge = badges.get('Teacher')
    func = models.Post.objects.filter(author=user, score__gt=1).count
    if check(badge, awards, func=func):
        create(request, user=user, badge=badge)
        return
    
    badge = badges.get('Supporter')
    func = models.Vote.objects.filter(author=user).count
    if check(badge, awards, func=func):
        create(request, user=user, badge=badge)
        return
    
    badge = badges.get('Nice Question')
    func = models.Post.objects.filter(author=user, score__gt=10).count
    if check(badge, awards, func=func):
        create(request, user=user, badge=badge)
        return
