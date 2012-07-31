from django.conf import settings
from main.server import models, html
from main.server.const import *

from django.db.models import Avg, Max, Min, Count

def remove_notes(target, maxcount=1000):
    """Clears the notes  for each user"""
    count = models.Note.objects.filter(target=target).count()
    if count > maxcount:   
        last  = models.Note.objects.filter(target=target).order_by('-date')[maxcount]
        clear = models.Note.objects.filter(target=target, date__lt=last.date).exclude(sender=target)
        clear.delete()
        count = models.Note.objects.filter(target=target).count()
        print '*** deleting for user %s, count %s' % (target.id, count)
    
def trim_notelist(maxcount=1000):
    query = models.User.objects.annotate(note_count=Count('note_target')).filter(note_count__gt = 2 * maxcount)
    for user in query:
        remove_notes(user, maxcount=maxcount)
    
def run():
    trim_notelist(maxcount=1000)
    
if __name__ == '__main__':
    run()