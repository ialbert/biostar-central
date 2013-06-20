"""
Bumps semi-random posts
"""
import random, datetime
    
from django.conf import settings
from main.server import models, html
from main.server.const import *

def bump():

    query = models.Post.objects.filter(type=POST_QUESTION, status=POST_OPEN)

    if random.random() > 0.75:
        query = query.filter(answer_count=0)

    query = query.values_list("id")

    ids = [ p[0] for p in query ]

    pk = random.choice(ids)

    community = models.User.objects.get(pk=1)
    post = models.Post.objects.get(pk=pk)
    post.lastedit_date = datetime.datetime.now()
    post.lastedit_user = community
    post.save()

    print '*** bumped %s' % post.title


if __name__ == '__main__':
    bump()