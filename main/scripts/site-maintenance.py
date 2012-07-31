from django.conf import settings
from main.server import models, html
from main.server.const import *

from django.db.models import Avg, Max, Min, Count

def clear_notes(target, maxcount=1000):
    """Clears the notes  for each user"""
    count = models.Note.objects.filter(target=target).count()
    if count > maxcount:   
        last  = models.Note.objects.filter(target=target).order_by('-date')[maxcount]
        clear = models.Note.objects.filter(target=target, date__lt=last.date)
        clear.delete()
        count = models.Note.objects.filter(target=target).count()
        print '*** deleting for user %s, count %s' % (target.id, count)
    
def run(maxcount=1000):
    query = models.User.objects.annotate(note_count=Count('note_target')).filter(note_count__gt = 2 * maxcount)
    for user in query:
        clear_notes(user, maxcount=maxcount)
    
if __name__ == '__main__':
    run(10)