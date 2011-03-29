from django import template
from biostar.server import models

register = template.Library()

def votebox(context, post):
    return { 'post':post,
            'upvoted': post.get_vote(context['user'], models.VOTE_UP) is not None,
            'downvoted':post.get_vote(context['user'], models.VOTE_DOWN) is not None}



register.inclusion_tag('widgets/votebox.html', takes_context=True)(votebox)





