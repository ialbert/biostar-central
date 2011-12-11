"""
Indexes all post content
"""    
    
import shutil
from whoosh import index
from django.conf import settings
from main.server import  models
from main.server.const import *
from itertools import *

if __name__ == '__main__':
    shutil.rmtree(settings.WHOOSH_INDEX)
    models.create_index()
    ix = index.create_in(settings.WHOOSH_INDEX, models.WhooshSchema)
    wr = ix.writer()

    print "*** whoosh indexing %s posts" % models.Post.objects.all().count()
    
    for step, post in izip(count(1), models.Post.objects.all()):
        if post.type in POST_CONTENT_ONLY:
            text = post.content
        else:
            text = post.title + post.content
        wr.add_document(content=text, pid=post.id)

        if step % 1000 == 0:
            wr.commit()
            wr = ix.writer()
    wr.commit()