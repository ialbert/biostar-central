"""
Indexes all post content
"""    
    
import shutil
from whoosh import index
from django.conf import settings
from main.server import  models, const

if __name__ == '__main__':
    shutil.rmtree(settings.WHOOSH_INDEX)
    models.create_index()
    ix = index.create_in(settings.WHOOSH_INDEX, models.WhooshSchema)
    wr = ix.writer()

    print "*** whoosh indexing %s posts" % models.Post.objects.all().count()
    for post in models.Post.objects.all():
        if post.post_type in const.POST_FULL_FORM:
            text = post.title + post.content
        else:
            text = post.content
        wr.add_document(content=text, pid=post.id)
    wr.commit()