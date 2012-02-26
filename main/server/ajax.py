#
# handler methods that handle all ajax based interactions
#
from functools import partial
from collections import defaultdict
from main.server import html, models, notegen, auth
from main.server.html import get_page
from main.server.const import *
from django.contrib.auth.decorators import login_required
from datetime import datetime, timedelta
from django.http import HttpResponse

# activate logging
import logging
logger = logging.getLogger(__name__)

def ajax_msg(msg, status):
    return html.json_response(dict(status=status, msg=msg))
    
ajax_success = partial(ajax_msg, status='success')
ajax_error   = partial(ajax_msg, status='error')

class ajax_error_wrapper(object):
    "used as decorator to trap/display  errors in the ajax calls"
    def __init__(self, f):
        self.f = f
        
    def __call__(self, *args, **kwds):
        try:
            # first parameter is the request
            if args[0].method != 'POST':
                return ajax_error('POST method must be used.')
               
            if not args[0].user.is_authenticated():
                return ajax_error('User must be logged in.')
                
            value = self.f(*args, **kwds)
            return value
        except Exception,exc:
            return ajax_error('Error: %s' % exc)

@ajax_error_wrapper           
def vote(request):
    "Handles all voting on posts"
    
    author = request.user
    
    # attempt to find the post and vote
    post_id = int(request.POST.get('post'))
    post = models.Post.objects.get(id=post_id)

    # get vote type
    type = request.POST.get('type')
    
    # remap to actual type
    type = dict(upvote=VOTE_UP, accept=VOTE_ACCEPT, bookmark=VOTE_BOOKMARK).get(type)
        
    if not type:
        return ajax_error('invalid vote type')
            
    if type  == VOTE_UP and post.author == author:
        return ajax_error('You may not vote up on your own post')
    
    if type == VOTE_ACCEPT and post.author == author:
        return ajax_error('You may not accept your own post')
    
    if type == VOTE_ACCEPT and post.root.author != author:
        return ajax_error('Only the original poster may accept an answer')
        
    # see if there is an existing vote of this type
    old_vote = post.get_vote(author, type)

    if old_vote:
        msg = '%s removed' % old_vote[0].get_type_display()
        post.remove_vote(author, type)
        logger.info('%s\t%s\t%s' % (author.id, post.id, msg) )
        return ajax_success(msg)
    
    if type == VOTE_BOOKMARK:
        vote = post.add_vote(author, type)
        return ajax_success('%s added' % vote.get_type_display())
    
    # throttle
    today  = datetime.now()
    shift  = timedelta(seconds=VOTE_SESSION_LENGTH)
    past   = today - shift
    count  = models.Vote.objects.filter(author=author, date__gt=past).count()
    avail  = MAX_VOTES_PER_SESSION - count
    
    if avail <= 0:
        msg = "You ran out of votes ;-) there will more in a little while"
        logger.info('%s\t%s\t%s' % (author.id, post.id, "out of votes") )
        return ajax_error(msg)
    else:
        # log all voting into the server log
        vote = post.add_vote(author, type)
        msg  = '%s added. %d votes left' % (vote.get_type_display(), avail-1)
        logger.info('%s\t%s\t%s' % (author.id, post.id, msg) )
        return ajax_success(msg)

@ajax_error_wrapper 
def moderate_post(request, pid, action):
    
    user = request.user
    
    action_map = { 'close':REV_CLOSE, 'reopen':REV_REOPEN,
        'delete':REV_DELETE, 'undelete':REV_UNDELETE }
    
    action_val = action_map.get(action)
    if not action_val:
        return ajax_error('Unrecognized action')
    
    post = models.Post.objects.get(id=pid)
    
    # two conditions where a post moderation may occur
    cond1  = user.can_moderate
    cond2  = (user == post.author) and action_val in (REV_CLOSE, REV_DELETE)
    permit = cond1 or cond2
    
    if permit :
        models.moderate_post(post=post, action=action_val, user=user)
        msg = '%s performed' % action
        return ajax_success(msg)
    
    return ajax_error('permission denied - please ask a moderator')
        
@ajax_error_wrapper 
def moderate_user(request, uid, action):
 
    user = request.user
    target = models.User.objects.get(id=uid)

    if not auth.authorize_user_edit(target=target, user=user):
        return ajax_error('Permission denied')
    
    if (target == user):
        return ajax_error('Users may not moderate themselves')
        
    action_map = { 'suspend':USER_SUSPENDED, 'reinstate':USER_ACTIVE }
    action_val = action_map.get(action)
    if not action_val:
        return ajax_error('Unrecognized action')
    
    msg = models.moderate_user(user=user, target=target, action=action_val)
    msg = "%s" % action.title()
    return ajax_success(msg)

@ajax_error_wrapper           
def post_destroy(request, pid):
    
    
    user = request.user
    post = models.Post.objects.get(id=pid)
    
    # two conditions where a comment destruction may occur
    permit  = user.can_moderate or (user == post.author )
    if not permit:
        return ajax_error('Permission denied')
        
    status = models.destroy_post(post=post, user=user)
    if status:
        return ajax_success('deleted')
    else:
        return ajax_success('destroyed')
    
    
@ajax_error_wrapper
def preview(request):
    "This runs the markdown preview functionality"
    content = request.POST.get('content','no input')[:5000]

    try:
        output = html.generate(content)
    except KeyError, exc:
        # return more userfriendly errors, used for debugging
        output = 'Error: %s' % str(exc)

    return HttpResponse(output, mimetype='text/plain')