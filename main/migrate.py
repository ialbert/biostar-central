"""
Parses a SE biostar xml datadump, and populates the 
database with the content. Finally exports a loadable 
data fixture from this database.
"""
import sys, os, random, re, shutil, gc, string
from datetime import datetime
from itertools import *
from xml.etree import ElementTree
from collections import defaultdict

#  we will need to fixup path so we can access the Django settings module
join = os.path.join
FILE_DIR = os.path.dirname(__file__)

# add paths and default settings module if nothing was set
if not os.environ.get('DJANGO_SETTINGS_MODULE'):
    sys.path.append( join(FILE_DIR, '..' ))
    sys.path.append( join(FILE_DIR, '..', 'main'))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

# load up the django framework with settings
from django.conf import settings
from main.server import models, const
from django.db import transaction
from django.db.models import signals
from django.utils.datastructures import SortedDict

# import all constants
from main.server.const import *

def xml_reader(fname, limit=None):
    """
    SE XML dumps use similar format, everything is in tags, no attributes used anywhere
    This iterator will parse the XML and return a list of dictionaries keyed by tag name and
    having the node value as the value of the key.
    
    Example: row><a>1</a><b>2</b></row> --> [ { 'a':'1', 'b':'2' } ]
    
    """
    elems = ElementTree.parse(fname)
    elems = islice(elems.findall('row'), limit)

    def keyfunc(x):
        return x.tag, x.text

    rows = []
    # each element is an iterable containing tags
    for elem in elems:
        tagmap = dict(map(keyfunc, elem))
        rows.append(tagmap)
    
    print '*** parsed %d records from %s' % (len(rows), fname)
    return rows

def parse_time(timestr):
    "Attempts to convert a timestring to a datetime object"
    try:
        return datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        return datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S')
    
# matches tags in a tag string
tag_finder = re.compile('[^a-z]([a-z]+)[^a-z]')
        
# returns a list of tags from a tagstring
def parse_tag_string(rawtagstr):
    """
    Parses tags in a SE dump
    """
    tags = [ t.strip() for t in tag_finder.findall(rawtagstr) ]
    tags = filter(None, tags)
    return ' '.join(tags)

def checkuser(row, key='Id'):
    """
    Checks a user for validity
    """
    userid = row.get(key)
    if not userid:
        if LOGSTREAM:
            out = [ '* attribute missing', key, row ]
            print '\t'.join(map(str, out))
        return False
    
    if userid == '-1':
        if LOGSTREAM:
            out = [ '* invalid user', key, row ]
            print '\t'.join(map(str, out))
        return False

    return True

@transaction.commit_manually
def insert_users(fname, limit):
    "Inserts the users"
    gc.collect()

    # needs to create community user
    user, flag = models.User.objects.get_or_create(username='community')
    user.profile.display_name = 'Biostar'
    user.profile.save()
    
    # reads the whole file
    rows = xml_reader(fname, limit=limit)
    
    rows = filter(checkuser, rows)

    # StackExchange to BioStar mappings
    typemap = { '4':const.USER_MODERATOR, '5':const.USER_ADMIN }

    # stores users and profiles
    ulist, plist = [], []

    # attribute collector
    for row in rows:
    
        # these are user related attributes
        userid   = row['Id']
        
        # de
        if userid == '-1':
            continue
        username = 'u%s' % userid
        email    = row.get('Email') or username
        name     = row.get('DisplayName', 'User %s' % userid)
        last_login  = parse_time(row['LastAccessDate'])
        date_joined = parse_time(row['CreationDate'])

        # stores the users
        upair = (userid, dict( username=username, first_name=name, 
                last_name="", email=email, last_login=last_login, date_joined=date_joined))
        ulist.append( upair ) 

        # these will be profile related attributes
        score = int(row.get('Reputation', 0))
        utype = typemap.get( row['UserTypeId'], const.USER_MEMBER)
        display_name = row.get('DisplayNameCleaned', 'User %s' % userid).title()
        website  = row.get('WebsiteUrl', '')
        about_me = row.get('AboutMe', '')
        location = row.get('Location', '')
        openid   = row.get('OpenId', 'http://www.biostars.org')
        last_login_ip = row.get('LastLoginIP', '0.0.0.0')

        # store profiles
        ppair = (userid, dict( display_name=display_name, website=website, location=location, 
                last_login_ip=last_login_ip, about_me=about_me, type=utype, last_visited=last_login ))        
        plist.append(ppair)

    users = {}
    print "*** inserting %s users" % len(ulist)
    with transaction.commit_on_success():
        for (userid, u), (userid2, p) in zip(ulist, plist):
            assert userid == userid2, 'Sanity check'
            user = models.User(**u) 
            if USE_DB:
                user.is_staff = p['type'] == const.USER_ADMIN
                user.is_superuser = p['type'] == const.USER_ADMIN                
                user.save()
                prof = user.get_profile()
                for attr, value in p.items():
                    setattr(prof, attr, value)
                prof.save()

              
            users[userid] = user

    return users

def checkfunc(key, data):
    """
    Returns a filtering fuction that safely checks 
    that the keys of their parameters are in an data structure 
    """
    global LOGSTREAM

    def func(row):
        cond1 = (key in row)
        if not cond1:
            if LOGSTREAM:
                out = [ '* attribute error', key, row ]
                print '\t'.join(map(str, out))
            return False

        cond2 = (row[key] in data)
        if not cond2:
            if LOGSTREAM:
                out = [ '* lookup error', key, row[key], row ]
                print '\t'.join(map(str, out))

            return False
        return True

    return func

def insert_posts(fname, limit, users):
    "Inserts the posts"
    gc.collect()

    # read all the posts
    rows = xml_reader(fname, limit=limit)

    # keep only posts with a valid user
    rows = filter(checkfunc('OwnerUserId', users), rows)

    plist   = [] # collect post attributes
    acount  = defaultdict(int) # maintains answer counts
    parents = dict() # maps questions to answers
    # first insert all posts
    
    # connects the post type in the SE dump to the models in BioStar
    PMAP = { '1': const.POST_QUESTION, '2': const.POST_ANSWER }
    
    for row in rows:
        postid = row['Id']
        views  = row['ViewCount']
        creation_date = parse_time(row['CreationDate'])
        author   = users[row['OwnerUserId']]
        parentid = row.get('ParentId')
        title  = row.get('Title','')
        tag_string = row.get('Tags', '')
        if tag_string:
            tag_string = parse_tag_string(tag_string)

        ptypeid   = row['PostTypeId']
        post_type = PMAP.get(ptypeid, const.POST_OTHER)
        
        # collect answercounts
        if post_type == POST_ANSWER and parentid:
            acount [parentid] += 1
            parents[postid] = parentid
    
        ppair = (postid, dict(author=author, views=views, creation_date=creation_date, type=post_type, title=title, tag_val=tag_string))
        plist.append (ppair)
        
    posts = {}
    print "*** inserting %s posts" % len(plist)
   
    with transaction.commit_on_success():
        for i, (postid, p) in enumerate(plist):
            post = models.Post(**p)
            parentid = parents.get(postid)
            if USE_DB:
                if (i % 1000 == 0):
                    print "*** commit at %s" % i
                    transaction.commit()
                # gets triggered as a signal
                #post.answer_count = acount.get(postid, 0)
                post.answer_count = 0
                if parentid:
                    parent = posts.get(parents[postid])
                    if not parent:
                        continue
                    post.parent = post.root = parent
                post.save()
                post.set_tags()
                
            posts[postid] = post
            
    return posts

@transaction.commit_manually
def insert_post_revisions(fname, limit, users, posts):
    """
    Inserts post revisions. Also responsible for parsing out closed/deleted 
    states from the post history log
    """
    gc.collect()

    # no limits are necessary since it is limited by posts and users already
    rows = xml_reader(fname, limit=None)
    revs = {} # Dictionary for fast GUID lookup
    ords = [] # Ensures order doesn't get mixed up
    
    # Stack overflow splits up modifications to title, tags, and content
    # as separate rows in the post history XML distinguished by type
    # We need to first go through and collect them together
    # based on the GUID they set on them.
       
    # keep revisions to valid posts/users
    rows = filter(checkfunc('UserId', users), rows)
    rows = filter(checkfunc('PostId', posts), rows)
    
    revs  = {}
    glist = [] # maintains order between revisions
    alist = [] # action list

    for row in rows:
        guid = row['RevisionGUID']
        date = parse_time(row['CreationDate'])
        post = posts[ row['PostId'] ]
        author = users[ row['UserId'] ]
        rtype  = row['PostHistoryTypeId']
        if guid not in revs:
            glist.append(guid)

        # we will update a revision to contain all changes
        rev = revs.get(guid, {'post':post, 'lastedit_user':author, 'lastedit_date':date})

        if rtype in ['1', '4']: # Title modification
            rev['title'] = row['Text']
        if rtype in ['2', '5']: # Body modification
            rev['content'] = row['Text']                     
        if rtype in ['3', '6']: # Tag modification
            rev['tag_string'] = parse_tag_string(row['Text'])

        if rtype in ['10','11','12','13']: # Moderator actions
            actions = {'10':const.POST_CLOSED, '11':const.POST_OPEN,
                       '12':const.POST_DELETED, '13':const.POST_OPEN}
            # this is defined in the models
            alist.append( (post, actions[rtype], author, date) )

        revs[guid] = rev

    print "*** inserting %s revisions" % len(glist)
    with transaction.commit_on_success():
        for (i, guid) in enumerate(glist):
            data = revs[guid]
            post = data['post']
            del data['post']
            if USE_DB:
                if (i % 1000 == 0):
                    print "*** commit at %s" % i
                    transaction.commit()
                for key, value in data.items():
                    setattr(post, key, value)
                post.save()

    print "*** inserting %s moderator actions" % len(alist)
    with transaction.commit_on_success():
        for post, status, user, date in alist:
            if USE_DB and post.id:
                models.post_moderate(post=post, status=status, user=user, date=date)
                
    # some posts may have been removed
    for (key, post) in posts.items():
        if not post.id:
            del posts[key]
    
def insert_votes(fname, limit, users, posts):

    gc.collect()
    rows = xml_reader(fname)

    rows = filter(checkfunc('UserId', users), rows)
    rows = filter(checkfunc('PostId', posts), rows)
    
    vlist = []
    for row in rows:
        post = posts[row['PostId']]
        user = users[row['UserId']]
        addr = users.get(row['IPAddress'])
        VoteType = row['VoteTypeId']
            
        if VoteType == '1':
            vote_type = const.VOTE_ACCEPT
        elif VoteType == '2':
            vote_type = const.VOTE_UP
        elif VoteType == '3':
            vote_type = const.VOTE_DOWN
        else:
            continue
        param = dict(post=post, author=user, type=vote_type)
        vlist.append(param)

    print "*** inserting %s votes" % len(vlist)
    with transaction.commit_on_success():
        for i, param in enumerate(vlist):
            vote = models.Vote(**param)
            if USE_DB:
                if (i % 1000 == 0):
                    print "*** commit at %s" % i
                    transaction.commit()
                vote.save()
        
def insert_comments(fname, posts, users, limit):
    
    gc.collect()

    rows = xml_reader(fname)

    # keep the valid rows only
    rows = filter(checkfunc('UserId', users), rows)
    rows = filter(checkfunc('PostId', posts), rows)
    
    comms, clist = {}, []
    for row in rows:
        cid  = row['Id']
        text = row['Text']
        postid = row['PostId'] 
        author = users[row['UserId']]
        creation_date = parse_time(row['CreationDate'])
        post_type = const.POST_COMMENT
        row = postid, cid, dict(author=author, creation_date=creation_date, content=text, type=post_type)
        clist.append( row )

    print "*** inserting %s comments" % len(clist)
    with transaction.commit_on_success():   
        for i, (postid, cid, param) in enumerate(clist):   
            parent = posts[postid]
            param['parent'] = parent
            param['root']   = parent.root or parent
            post = models.Post(**param)
            comms[cid] = post
            if USE_DB:
                if (i % 1000 == 0):
                    print "*** commit at %s" % i
                    transaction.commit()
                    gc.collect()
                post.save() 

    return comms

def insert_comment_votes(fname, limit, comms, users):
    "Inserts vote on comments"
    
    gc.collect()

    rows = xml_reader(fname)
    rows = filter(checkfunc('UserId', users), rows)
    rows = filter(checkfunc('PostCommentId', comms), rows)
    
    vlist = []
    for row in rows:
        post = comms[row['PostCommentId']]
        user = users[row['UserId']]
        addr = users.get(row['IPAddress'])
        VoteType = row['VoteTypeId']
            
        if VoteType == '1':
            vote_type = const.VOTE_ACCEPT
        elif VoteType == '2':
            vote_type = const.VOTE_UP
        elif VoteType == '3':
            vote_type = const.VOTE_DOWN
            # we will not migrate over downvotes
            continue
        else:
            continue
        if not post.id:
            continue
        param = dict(post=post, author=user, type=vote_type)
        vlist.append(param)

    print "*** inserting %s comment votes" % len(vlist)
    with transaction.commit_on_success():
        for (i, param) in enumerate(vlist):
            vote = models.Vote(**param)
            if USE_DB:
                if (i % 1000 == 0):
                    print "*** commit at %s" % i
                    transaction.commit()
                vote.save()
 
def insert_badges(fname, limit):
    "Inserts the badges"

    gc.collect()

    blist = []
    rows = xml_reader(fname, limit=limit)
    bmap = {}
    BCONST = { '3':models.BADGE_BRONZE, '2':models.BADGE_SILVER, '1':models.BADGE_GOLD }

    for row in rows:
        bid    = row['Id']
        bclass = row['Class']
        name   = row['Name']
        desc   = row['Description']
        unique = row['Single'] == 'true'
        secret = row['Secret'] == 'true'
        btype  = BCONST[bclass]
        param = dict(name=name, type=btype, description=desc, unique=unique, secret=secret)
        blist.append((bid, param))
            
    print "*** inserting %s badges" % len(blist)
    with transaction.commit_on_success():   
        for (bid, param) in blist:
            badge = models.Badge(**param)
            if USE_DB:
                badge.save()

            bmap[bid] = badge    
    return bmap

def finalize():
    print "*** finalizing"
    users = models.User.objects.filter().all()
    community = models.User.objects.get(username='community')
    with transaction.commit_on_success():   
        # create a welcome message
        for user in users:
            models.Note.objects.create(sender=community, target=user, content='Welcome to **Biostar!**', type=const.NOTE_USER)

def insert_awards(fname, users, badges, limit):
    "Inserts the badge awards"
    
    rows = xml_reader(fname)
    rows = filter(checkfunc('UserId', users), rows)
    rows = filter(checkfunc('BadgeId', badges), rows)

    alist = []
    for row in rows:
        Id = row['Id']
        date = parse_time(row['Date'])
        user = users[row['UserId']]
        badge = badges[row['BadgeId']]
        param = dict(user=user, badge=badge, date=date)
        alist.append(param)
    
    print "*** inserting %s awards" % len(alist)
    with transaction.commit_on_success(): 
        for param in alist:
            award = models.Award(**param)
            if USE_DB:
                award.save()
    if USE_DB:
        models.Badge.objects.create(name=const.BETA_TESTER_BADGE , 
            description="Participated in the BioStar public beta test. Thank you!", type=const.BADGE_GOLD, unique=True)

def admin_init():
    
    # create the admin users for a given settings file
    for index, (name, email) in enumerate(settings.ADMINS):
        
        admins = models.User.objects.filter(email=email).all()
            
        if not admins:
            # create this adminsitrative user
            admin  = models.User(email=email, first_name=name, username='a%i' % index)
            admins = [ admin ]
            print '*** created admin user %s (%s)' % (admin.username, admin.email)

        for admin in admins:
            # add admin rights to each user
            admin.set_password(settings.SECRET_KEY)
            admin.is_staff = admin.is_superuser = True
            if USE_DB:
                admin.save()
                admin.profile.type = const.USER_ADMIN
                admin.profile.save()
            print '*** added staff access to admin user %s (%s)' % (admin.username, admin.email)
    
    editors, flag = models.Group.objects.get_or_create(name=const.MODERATOR_GROUP)
    if flag:
        print '*** created group %s' % editors.name


def parse_post(fname):
    "An importable post is self contained with title and tags in it"
    lines = file(fname).readlines()
    title = lines[0].split(':')[-1]
    tags  = lines[1].split(':')[-1]
    body  = ''.join(lines[2:])
    return map(string.strip, (title, tags, body))

def tutorial_init():
    from glob import glob

    # find the first admin user
    admin = models.User.objects.get(email=settings.ADMINS[0][1])
    

    for fname in glob('import/forum/*.txt'):
        print "*** importing forum post %s" % fname
        title, tag_val, body = parse_post(fname)
        post = models.Post(title=title, author=admin,  type=POST_FORUM, tag_val=tag_val, content=body)
        post.save()
        post.set_tags()
    
    user, flag = models.User.objects.get_or_create(username='angus')
    user.profile.display_name = "MSU course 2011"
    user.profile.website  = "http://ged.msu.edu/angus/tutorials-2011/"
    user.profile.about_me = 'ANGUS is a site built around the `2010 course on [Analyzing Next-Generation Sequencing Data](http://bioinformatics.msu.edu/ngs-summer-course-2010 "Angus Website")'
    user.profile.save()
    for fname in glob('import/tutorials/2011/*.rst'):
        print "*** importing tutorial post %s" % fname
        title, tag_val, body = parse_post(fname)
        post = models.Post(title=title, author=user,  type=POST_TUTORIAL, tag_val=tag_val, content=body)
        post.save()
        post.set_tags()
        
def execute(path, limit=None):
    """
    Executes the imports
    """
    
    # insert users into the database
    fname = join(path, 'Users.xml')
    
    # no real user file, use anonymized users
    if not os.path.isfile(fname):
        print '*** using the anonymized user file'
        fname = join(path, 'AnonUsers.xml')
        
    users = insert_users(fname=fname, limit=limit)
  
    #tutorial_init()
    #return

    fname = join(path, 'Posts.xml')
    posts = insert_posts(fname=fname, limit=limit, users=users)
   
    # this creates way too many notices, so disconnect it during imports
    signals.post_save.disconnect( models.post_create_notification, sender=models.Post )

    fname = join(path, 'PostHistory.xml')
    revisions = insert_post_revisions(fname=fname, limit=limit, posts=posts, users=users)
       
    fname = join(path, 'PostComments.xml')
    comms = insert_comments(fname=fname, posts=posts, users=users, limit=limit)
    
    fname = join(path, 'Posts2Votes.xml')
    insert_votes(fname=fname, limit=limit, posts=posts, users=users)
    
    fname = join(path, 'Comments2Votes.xml')
    insert_comment_votes(fname=fname, limit=limit, comms=comms, users=users)
    
    fname = join(path, 'Badges.xml')
    badges = insert_badges(fname=fname, limit=limit)
    
    fname = join(path, 'Users2Badges.xml')
    insert_awards(fname=fname, users=users, badges=badges, limit=limit)
    
    # adds administration rights to users
    # listed in the DJAGNO settings file
    admin_init()
    tutorial_init()
    finalize()

if __name__ =='__main__':
    
    import doctest, optparse
    global LOGSTREAM, USE_DB

    # for debugging
    #sys.argv.extend( ["-p", "se0"] )
    
    # options for the program
    parser = optparse.OptionParser()
    parser.add_option("-o", dest="output", help="outputfile", default="import/datadump.json")
    parser.add_option("--path", dest="path", help="directory or zip archive containing a full biostar SE1 datadump")
    parser.add_option("--limit", dest="limit", help="limit to these many rows per file", type=int, default=None)
    parser.add_option("--log", dest="log", help="print information on entries that the migration skips", action="store_true", default=False)
    parser.add_option("--dry", dest="dry", help="dry run, no database action, checks the formats", action="store_true", default=False)
   
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not opts.path:
        parser.print_help()
        sys.exit()
        
    # also run the doctests
    #doctest.testmod()
    
    print '*** migration path %s' % opts.path
    if opts.log:
        LOGSTREAM = sys.stdout 
    else:
        LOGSTREAM = None
    
    # wether to use the database
    USE_DB = not opts.dry

    # call into the main program
    from django.core import management
    management.call_command('syncdb', verbosity=0, interactive=False)

    execute(path=opts.path, limit=opts.limit)
    
    fp = file(opts.output, 'wt')
    sys.stdout = fp
    management.call_command('dumpdata', 'auth.User', 'server', verbosity=1, interactive=False)