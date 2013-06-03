"""
Biostar views
"""
import difflib, time, re, random
from datetime import datetime, timedelta

from functools import partial
from collections import defaultdict
from main.server import html, models, const, formdef, action, notegen, auth
from main.server.html import get_page
from datetime import datetime
from django.db import connection

from django.contrib.auth.decorators import login_required
from django.db import transaction
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.contrib.auth import authenticate, login, logout
from django.contrib import messages
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Q
# the openid association model

from django.core.urlresolvers import reverse
from django.db.models import Count, Max, Min

# import all constants
from main.server import const
from main.server.const import *
from main import middleware

# activate logging
import logging
logger = logging.getLogger(__name__)


def get_post_manager(request):
    user = request.user
    if user.is_authenticated() and user.profile.can_moderate:
        return models.Post.objects.select_related('author', 'author__profile')
    else:
        return models.Post.open_posts.select_related('author', 'author__profile')

# mapst a word to a numeric post type
POST_TYPE_MAP = dict(
    questions=POST_QUESTION,  tutorials=POST_TUTORIAL, answers=POST_ANSWER, videos=POST_VIDEO,
    planet=POST_BLOG, tools=POST_TOOL, jobs=POST_JOB, news=POST_NEWS, publications=POST_PUBLICATION,
    forum=POST_FORUM,
)

POST_TYPE_REV_MAP = dict( [ (v,k) for (k,v) in POST_TYPE_MAP.items()] )

def mytags_posts(request):
    "Gets posts that correspond to mytags settins or sets warnings"
    user = request.user
    if not user.is_authenticated():
        messages.warning(request, "This Tab is populated only for registered users based on the My Tags field in their user profile")
        text = ""
    else:
        text = request.user.profile.my_tags 
        if not text:
            messages.warning(request, "Showing posts matching the My Tags fields in your user profile. Currently this field is not set.")
            
    return models.query_by_tags(user, text=text)

def filter_by_date(request, posts, since):
    days = 0
    if since == "today":
        days = 1
    elif since == "this week":
        days = 7
    elif since == "this month":
        days = 30
    elif since == "this year":
        days = 365
    if days:
        now = datetime.now().date()
        start = now - timedelta(days=days)
        posts = posts.filter(lastedit_date__gt=start)

    return posts

def filter_by_type(request, posts, post_type):
    "Filters posts by type"

    # filter is a single type
    if post_type == POST_FORUM:
        return posts.filter(type__in=[POST_FORUM, POST_NEWS])
    elif post_type in POST_TYPE_REV_MAP:
        return posts.filter(type=post_type)
    elif post_type == "howto":
        return posts.filter(type__in=[POST_TUTORIAL, POST_TIP, POST_TOOL, POST_VIDEO, POST_REVIEW])
    elif post_type == "myvotes":
        return  posts.filter(votes__post__author=request.user)
    elif post_type == "myposts":
        return  posts.filter(author=request.user)
    elif post_type == "mybookmarks":
        return  posts.filter(votes__author=request.user, votes__type=const.VOTE_BOOKMARK)
    elif post_type == 'sticky':
        return posts.filter(sticky=True)
    elif post_type == 'unanswered':
        return posts.filter(type__in=[POST_QUESTION, POST_FIXME], answer_count=0)
    elif post_type == 'all':
        return posts.exclude(type__in=POST_EXCLUDE)
    elif post_type == 'mytags':
        return mytags_posts(request)
    elif post_type == 'recent':
        return posts.exclude(type=POST_BLOG).select_related('author', 'author__profile','root')
    elif post_type == "galaxy":
        return posts.filter(type__in=POST_TOPLEVEL, tag_set__name__in=["galaxy"])
    elif post_type == "bookmarked":
        return posts.filter(book_count__gt=0)
    return posts.exclude(type__in=POST_EXCLUDE)
    
def apply_sort(request, posts, order, sticky=True):
    "Sorts posts by an order"
    sort_order = SORT_MAP.get(order, "-lastedit_date")
    if sticky:
        args = [ "-sticky", sort_order]
    else:
        args = [ sort_order ]
    return posts.order_by(*args)

# there is a tab bar and a lower "pill" bar

SORT_CHOICES  = "activity,rank,views,votes,answers,edit".split(',')
DATE_FILTER   = "all time,today,this week,this month,this year".split(',')

def index(request, tab='all'):
    user = request.user
    auth = user.is_authenticated()
    
    # asking for an invalid tab
    if tab not in VALID_TABS:
        msg = html.sanitize('Unknown content type "%s"' % tab)
        messages.error(request, msg)
        return html.redirect("/")
        
    # populate the session data
    sess = middleware.Session(request)
    
    # get the sort order
    sort_type = sess.sort_order()

    # parse the date request
    since = request.GET.get('since', DATE_FILTER[0]).lower()

    # set the last active tab
    sess.set_tab(tab)
    
    # get the numerical value for these posts
    post_type = POST_TYPE_MAP.get(tab, tab)

    # override the sort order if the content so requires
    sort_type = 'creation' if tab=='recent' else sort_type

    # this here needs to be reworked TODO
    if tab == "best":
        sort_type = "votes"
        since = request.GET.get('since', 'this week')
        messages.info(request, "Most <b>upvoted</b> active posts of <b>%s!</b>" % since)

    elif tab == "bookmarked":
        sort_type = "bookmark"
        since = request.GET.get('since', 'this month')
        messages.info(request, "Most <b>bookmarked</b> active posts of <b>%s!</b>" % since)
    elif tab == "myvotes":
        sort_type = "votes__date"
        messages.info(request, "Posts created by you that have received up-votes from other users")
    elif tab == "mybookmarks":
        sort_type = "votes__date"
        messages.info(request, "Your bookmarked posts")
    elif tab == "myposts":
        sort_type = "creation"
        messages.info(request, "Posts created by you")

    # the params object will carry
    layout = settings.USER_PILL_BAR if auth else settings.ANON_PILL_BAR
    
    # wether to show the type of the post
    show_type = post_type in ('all', 'recent')

    show_search = True

    if tab in VALID_PILLS:
        tab, pill = "posts", tab
    else:
        tab, pill = tab, ""

    params  = html.Params(tab=tab, pill=pill, sort=sort_type, sort_choices=SORT_CHOICES, date_filter=DATE_FILTER, since=since,
                          layout=layout, show_type=show_type, title="Bioinformatics Answers", show_search=show_search)

    # this will fill in the query (q) and the match (m)parameters
    params.parse(request)
    
    # returns the object manager that contains all or only visible posts
    posts = get_post_manager(request)

    # filter posts by type
    posts = filter_by_type(request=request, posts=posts, post_type=post_type)

    # apply date filtering
    posts = filter_by_date(request=request, posts=posts, since=since)

    # reduce SQL query count by preselecting data that will be displayed
    posts = posts.select_related('author', 'author__profile', 'lastedit_user', 'lastedit_user__profile')
        
    # sticky is not active on recent and all pages
    sticky = (tab != 'recent') and (pill not in ('all', "best", "bookmarked"))

    # hackfix
    sticky = False if pill.startswith("my") else sticky

    # order may change if it is invalid search
    posts = apply_sort(request=request, posts=posts, order=sort_type, sticky=sticky)

    # get the counts for the session
    counts = sess.get_counts(post_type)

    page = get_page(request, posts, per_page=settings.POSTS_PER_PAGE)
    
    # save the session
    sess.save()

    # try to set a more informative title
    title_map = dict(
            questions="Bioinformatics Questions", unanswered="Unanswered Questions", tutorials="Bioinformatics Tutorials",
            jobs="Bioinformatics Jobs", videos="Bioinformatics Videos", news='Bioinformatics News', tools="Bioinformatics Tools",
            recent="Recent bioinformatics posts", planet="Bioinformatics Planet",
            galaxy="Galaxy on Biostar", bookmarked="Most bookmarked",
    )


    params.title = title_map.get(pill) or title_map.get(tab, params.title)

    return html.template(request, name='index.html', page=page, params=params, counts=counts)

# generate the Best Of tab, collection a highest rated posts
def bestof(request, tab='best'):
    messages.info(request, "Most <b>upvoted</b> active posts <b>this week!</b>")
    return html.redirect("/show/best/?sort=votes&since=this week")


def show_tag(request, tag_name=''):
    "Display posts by a certain tag"
    user = request.user
    # populate the session data
    sess = middleware.Session(request)
    
    # get the sort order
    sort_type = sess.sort_order()

    # parse the date request
    since = request.GET.get('since', DATE_FILTER[0]).lower()

    # select based on history
    tab, pill = "posts", sess.get_tab()
    
    params = html.Params(nav='', tab=tab, sort='' )
    
    # the params object will carry
    layout = settings.USER_PILL_BAR if auth else settings.ANON_PILL_BAR
    
    # wether to show the type of the post
    params  = html.Params(tab=tab, pill=pill, sort=sort_type, sort_choices=SORT_CHOICES, date_filter=DATE_FILTER, since=since,
            layout=layout, title="Tagged as %s" % tag_name)
    
    msg = 'Filtering by tag: <b>%s</b>. Subscribe to an <a href="/feeds/tag/%s/">RSS feed</a> to this tag.' % (tag_name,tag_name)
    messages.info(request, msg)
    posts = models.query_by_tags(user=user, text=tag_name)

    # apply date filtering
    posts = filter_by_date(request=request, posts=posts, since=since)

    # order may change if it is invalid search
    posts = apply_sort(request=request, posts=posts, order=sort_type, sticky=False)

    page  = get_page(request, posts, per_page=20)

    return html.template( request, name='index.html', page=page, params=params)

def show_user(request, uid, post_type=''):
    "Displays posts by a user"

    user = models.User.objects.filter(id=uid).select_related('profile').all()[0]
    params = html.Params(nav='', tab='user', sort='', title="Activity for user %s" % user.profile.display_name)

    # notification
    messages.info(request, 'Filtering by user: %s' % user.profile.display_name)
   
    post_type = POST_REV_MAP.get(post_type.lower())
    if post_type:
        posts = get_post_manager(request).filter(type=post_type, author=user).order_by('-creation_date')
    else:
        posts = get_post_manager(request).filter(author=user).order_by('-creation_date')

    posts = posts.select_related('author', 'author__profile', 'root')
    page  = get_page(request, posts, per_page=20)
    return html.template( request, name='index.html', page=page, params=params)

def user_profile_redirect(request, uid, tab='activity'):
    """
    User's profile page
    """
    url = reverse("user-profile", kwargs=dict(uid=uid))
    return html.redirect(url, permanent=True)

def user_profile(request, uid, tab='activity'):
    """
    User's profile page
    """

    if not models.User.objects.filter(id=uid):
        messages.error(request, "This user does not exist. It has perhaps been deleted.")
        return html.redirect("/")
        
    user = request.user
    target = models.User.objects.get(id=uid)
    awards = []
    page   = None
    
    # some information is only visible to the user
    target.writeable = auth.authorize_user_edit(target=target, user=user, strict=False)
    show_private = target.showall = (target == user)

    params = html.Params(tab=tab, sort='', title="User %s" % target.profile.display_name)

    # these do not actually get executed unless explicitly rendered in the page
    bookmarks = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).select_related('post', 'post__author__profile').order_by('id')
    awards = models.Award.objects.filter(user=target).select_related('badge').order_by('-date')

    # we need to collate and count the awards
    answer_count = models.Post.objects.filter(author=target, type=POST_ANSWER).count()
    question_count = models.Post.objects.filter(author=target, type=POST_QUESTION).count()
    comment_count = models.Post.objects.filter(author=target, type=POST_COMMENT).count()
    post_count = models.Post.objects.filter(author=target).count()
    vote_count = models.Vote.objects.filter(author=target).count()
    award_count = models.Award.objects.filter(user=target).count()
    note_count  = models.Note.objects.filter(target=target, unread=True).count()
    bookmarks_count  = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).count()

    note_types = [ NOTE_USER, NOTE_PRIVATE ] if show_private else [ NOTE_USER ]

    if tab in [ 'activity', 'created' ]:
        if tab == 'created':
            q = Q(sender=target, target=target, type__in=note_types)
            e = Q()
        else:
            q = Q(target=target, type__in=note_types)
            e = Q(sender=target, type=NOTE_USER)

        notes = models.Note.objects.filter(q)
        notes = notes.exclude(e)
        notes = notes.select_related('author', 'author__profile', 'root').order_by('-date')

        page  = get_page(request, notes, per_page=15)
        # we evalute it here so that subsequent status updates won't interfere
        page.object_list = list(page.object_list)
        if show_private:
            models.Note.objects.filter(target=target, unread=True).update(unread=False)
            models.UserProfile.objects.filter(user=target).update(new_messages=0)
            note_count = 0

    elif tab == 'bookmarks':
        bookmarks = models.Vote.objects.filter(author=target, type=VOTE_BOOKMARK).select_related('post', 'post__author__profile').order_by('-date')
        page  = get_page(request, bookmarks, per_page=15)

    elif tab =="moderator":
        notes = models.Note.objects.filter(target=target, type=NOTE_MODERATOR).select_related('author', 'author__profile').order_by('-date')
        page  = get_page(request, notes, per_page=15)

    elif tab =="votes":
        votes = models.Vote.objects.filter(post__author=target).select_related('author', 'author__profile', 'post').order_by('-date')[:30]
        page  = get_page(request, votes, per_page=15)

    elif tab =="supporters":
        votes = models.Vote.objects.filter(post__author=target).select_related('author', 'author__profile', 'post').order_by('-date')
        votes = list(votes)[:21]
        if len(votes) < 10:
            votes = []
        random.shuffle(votes)
        page  = get_page(request, votes, per_page=21)


    params.update(dict(question_count=question_count, answer_count=answer_count, note_count=note_count, bookmarks_count=bookmarks_count,
            comment_count=comment_count, post_count=post_count, vote_count=vote_count, award_count=award_count))
    
    return html.template(request, name='user.profile.html', awards=awards,
        user=request.user,target=target, params=params, page=page)

def user_list(request):
    search  = request.GET.get('m','')[:80] # trim for sanity
    params = html.Params(nav='users', sort='')
    if search:
        query = Q(profile__display_name__icontains=search)
    else:
        query = Q(id__gt=0)
        
    users = models.User.objects.filter(query).select_related('profile').order_by("-profile__score", "id")
    page  = get_page(request, users, per_page=28)
    return html.template(request, name='user.list.html', page=page, params=params)

def tag_list(request):
    
    # remove null tags
    models.Tag.objects.all().filter(count=0).delete()
    
    search  = request.GET.get('m','')[:80] # trim for sanity
    
    if search:
        query = Q(name__icontains=search)
    else:
        query = Q(id__gt=0)
        
    tags = models.Tag.objects.filter(query).exclude(name="deleted-post").order_by('name')
    page = get_page(request, tags, per_page=250)
    params = html.Params(nav='tags', sort='', since='')
    return html.template(request, name='tag.list.html', page=page, params=params)

def badge_list(request):
    user = request.user
    
    badges = models.Badge.objects.filter(secret=False).order_by('-count', '-type')
    
    # set a flag for badges that a user has
    if user.is_authenticated():
        earned = set( models.Award.objects.filter(user=user).values_list('badge__name', flat=True).distinct() )
    else:
        earned = []
    for badge in badges:
        badge.earned = badge.name in earned
        
    params = html.Params(nav='badges', sort='')
    return html.template(request, name='badge.list.html', badges=badges, params=params)

def post_show_redirect(request, pid):
    """
    Permanent redirect from an old style post show
    """
    url = reverse("post-show", kwargs=dict(pid=pid))
    return html.redirect(url, permanent=True)

def post_show(request, pid):
    """
    Displays a post and its children
    """
    user = request.user

    # populate the session data
    sess = middleware.Session(request)
    tab  = "posts"
    pill = sess.get_tab()
    
    auth = user.is_authenticated()
    layout = settings.USER_PILL_BAR if auth else settings.ANON_PILL_BAR
    
    params  = html.Params(tab=tab, pill=pill, layout=layout)
    
    query = get_post_manager(request)

    try:
        root = query.get(id=pid)
        if not root.top_level:
            return html.redirect( root.get_absolute_url() )

        # update the views for the question
        models.update_post_views(post=root, request=request, minutes=const.POST_VIEW_UPDATE)
        counts = sess.get_counts()
    
    except models.Post.DoesNotExist, exc:
        messages.warning(request, 'The post that you are looking for does not exists. Perhaps it was deleted!')
        return html.redirect("/")
    
    # get all answers to the root
    children = models.Post.objects.filter(root=root).exclude(type=POST_COMMENT, id=root.id).select_related('author', 'author__profile').order_by('-accepted', '-score', 'creation_date')
    
    # comments need to be displayed by creation date
    comments = models.Post.objects.filter(root=root, type=POST_COMMENT).select_related('author', 'author__profile').order_by('creation_date')

    all = [ root ] + list(children) + list(comments)
    # add the various decorators
    
    models.decorate_posts(all, user)
    
    # may this user accept answers on this root
    accept_flag = (user == root.author)
    
    # these are all the answers
    answers = [ o for o in children if o.type == POST_ANSWER ]
    for a in answers:
        a.accept_flag = accept_flag
        
    # get all the comments
    tree = defaultdict(list)
    for comment in comments:  
        tree[comment.parent_id].append(comment)
   
    # generate the tag cloud
    #tags = models.Tag.objects.all().order_by('-count')[:50]
    
    return html.template( request, name='post.show.html', root=root, answers=answers, tree=tree, params=params, counts=counts)
 
def redirect(post):
    return html.redirect( post.get_absolute_url() )
    
@login_required(redirect_field_name='/openid/login/')
def new_comment(request, pid=0):
    "Shortcut to new comments"
    return new_post(request=request, pid=pid, post_type=POST_COMMENT)

@login_required(redirect_field_name='/openid/login/')
def new_answer(request, pid):
    return new_post(request=request, pid=pid, post_type=POST_ANSWER)

def safe_context(key_name, data):
    try:
        patt = settings.EXTERNAL_AUTHENICATION[key_name][1]
        return patt % data
    except Exception, exc:
        return "context error: %s" % exc

@login_required(redirect_field_name='/openid/login/')
def new_post(request, pid=0, post_type=POST_QUESTION, data=None):
    "Handles the creation of a new post"
    
    user   = request.user
    name   = "post.edit.html"
    parent = models.Post.objects.get(pk=pid) if pid else None
    root   = parent.root if parent else None
    toplevel = (pid == 0)
    factory  = formdef.ChildContent if pid else formdef.TopLevelContent


    params = html.Params(tab='new', title="New post", toplevel=toplevel, data=data)

    extern = formdef.ExternalLogin(request.GET)
    if extern.is_valid():
        e = extern.cleaned_data['data']
        context = safe_context(extern.cleaned_data['name'], e)
        tag_val = e.get("tags", "")
        title   = e.get("title", "")
        form = factory(initial=dict(tag_val=tag_val, context=context, title=title))
        return html.template(request, name=name, form=form, params=params)

    if request.method == 'GET':
        # no incoming data, render form
        form = factory()
        return html.template(request, name=name, form=form, params=params)

    # polls can only have comments associated with the root
    is_poll = root and root.type == POST_POLL
    is_author = root and root.author == user
    is_root_comment = (post_type == POST_COMMENT) and (parent == root)

    if is_poll and not is_author and not is_root_comment:
        messages.error(request, "Polls are special content. Users other than the author of the poll may only create comments directly associated with the top post.")
        return redirect(root)

    # process the incoming data
    assert request.method == 'POST', "Method=%s" % request.method

    form = factory(request.POST)

    # applying some throttles here
    now = datetime.now()

    # minimum age
    trust_range = timedelta(hours=settings.TRUST_INTERVAL)

    # a user that joined six hours ago
    new_user = (now - user.date_joined) < trust_range

    # posts within a trust interval
    post_count = models.Post.objects.filter(author=user, creation_date__gt=now - trust_range).count()

    # how many posts are too many
    too_many = (post_count >= settings.TRUST_NEW_USER_MAX_POST) if new_user else post_count >= settings.TRUST_USER_MAX_POST

    # different error messages for different type of users
    if new_user and too_many:
        messages.error(request, """
            Brand new users (accounts less than %s hours old) may not create more than %s posts.<br>
            Apologies if that interferes with your usage, it is an anti-spam measure.<br>
            Gives us a chance to ban the spam bots before too long. This limitation is lifted later.
            """ % (settings.TRUST_INTERVAL, settings.TRUST_NEW_USER_MAX_POST))
        return html.template(request, name=name, form=form, params=params)

    # this is main user throttling
    if too_many:
        messages.error(request, """
            A user may only create %s posts within a six hour period.<br>
            This is to protect against runaway bots.<br>
            Sorry for any inconvenience this may cause.
            """ % settings.TRUST_USER_MAX_POST)
        return html.template(request, name=name, form=form, params=params)

    if not form.is_valid():
        # returns with an error message
        return html.template(request, name=name, form=form, params=params)

    # form is valid at this point, create the post
    params = dict(author=user, type=post_type, parent=parent, root=root)

    # form may contain a variable number of elements
    for attr in "title content tag_val type context".split():
        if attr in form.cleaned_data:
            params[attr] = form.cleaned_data[attr]

    with transaction.commit_on_success():
        post = models.Post.objects.create(**params)
        post.set_tags()

    return redirect(post)

@login_required(redirect_field_name='/openid/login/')
def post_edit(request, pid=0):
    "Handles the editing of an existing post"
    
    user = request.user
    name = "post.edit.html"
    post = models.Post.objects.get(pk=pid)

    if not post.open and not user.can_moderate:
        messages.error(request, 'Post is closed. It may not be edited.')
        return redirect(post.root)
        
    # verify that this user may indeed modify the post
    if not auth.authorize_post_edit(post=post, user=request.user, strict=False):
        messages.error(request, 'User may not edit the post.')
        return redirect(post.root)
  
    toplevel = post.top_level
    factory  = formdef.TopLevelContent if toplevel else formdef.ChildContent

    params = html.Params(tab='edit', title="Edit post", toplevel=toplevel)
    if request.method == 'GET':
        # no incoming data, render prefilled form
        form = factory(initial=dict(title=post.title, content=post.content, tag_val=post.tag_val, type=post.type, context=post.context))
        return html.template(request, name=name, form=form, params=params)

    # process the incoming data
    assert request.method == 'POST', "Method=%s" % request.method
    form = factory(request.POST)
    if not form.is_valid():
        # returns with an error message
        return html.template(request, name=name, form=form, params=params)

    # form is valid now set the attributes
    for key, value in form.cleaned_data.items():
        setattr(post, key, value)
        post.lastedit_user = user
        post.lastedit_date = datetime.now()
        post.set_tags() # this saves the post
    
    return redirect(post)
    
def revision_show(request, pid):
    post = models.Post.objects.get(pk=pid)
    revs = post.revisions.order_by('-date').select_related('author', 'lastedit_author')
    return html.template(request, name='revision.show.html', revs=revs, post=post)
   
def post_redirect(request, pid):
    "Redirect to a post"
    post = models.Post.objects.get(id=pid)
    return html.redirect( post.get_absolute_url() )

def linkout(request, pid):
    "Used to be able to count the views for a linkout"
    post = models.Post.objects.get(id=pid)
    models.update_post_views(post=post, request=request)
    post.set_rank()
    post.save()
    if post.url:
        return html.redirect(post.url)    
    else:
        messages.error(request, 'linkout used on a post with no url set %s' % post.id)
        return html.redirect("/")

def modlog_list(request):
    "Lists of all moderator actions"
    mods = models.Note.objects.filter(type=NOTE_MODERATOR).select_related('sender', 'target', 'post', 'sender_profile').order_by('-date')
    page = get_page(request, mods)
    return html.template(request, name='modlog.list.html', page=page)

