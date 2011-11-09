"""
Parses a SE biostar xml datadump, and populates the 
database with the content. Finally exports a loadable 
data fixture from this database.
"""
import sys, os, random, re
from datetime import datetime
from itertools import *
from xml.etree import ElementTree

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
from main.server import models, const, siteinit
from django.db import transaction

def xml_reader(fname, limit=None):
    """
    SE XML dumps use similar format, everything is in tags, no attributes used anywhere
    This iterator will parse the XML and return a list of dictionaries keyed by tag name and
    having the node value as the value of the key.
    
    Example: row><a>1</a><b>2</b></row> --> [ { 'a':'1', 'b':'2' } ]
    
    """
    elems = ElementTree.parse(fname)
    elems = islice(elems.findall('row'), limit)

    rows = []
    for elem in elems:
        # transforms children nodes to a dictionary keyed by tag with the node text as value"
        pairs = map( lambda x: (x.tag, x.text), elem)
        data  = dict(pairs)
        rows.append( data )

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
    # parse_tag_string(" xtagx ybagy zmagz ")
    # 'tag1 tag2 tag3
    """
    tags = [ t.strip() for t in tag_finder.findall(rawtagstr) ]
    tags = filter(None, tags)
    return ' '.join(tags)

@transaction.commit_manually
def insert_users(fname, limit):
    "Inserts the users"
    
    store = {}
    rows = xml_reader(fname, limit=limit)
    
    # two step process generate the users then update the profiles
    for (index, row) in enumerate(rows):
        
        userid   = row['Id']
        username = 'u%s' % userid
        email    = row.get('Email') or username
        name     = row.get('DisplayName', username)
        try:
            user, flag = models.User.objects.get_or_create(username=username, email=email, first_name=name)
        except Exception, e:
            print 'Failed inserting row %s' % row
            print userid, name, username
            raise(e)
        
        # original ids will connect users to posts
        store[ userid ] = (flag, user, row)
    
    # these are the new users
    newusers = filter(lambda row: row[0], store.values())
    
    # this commit the user inserts
    print "*** inserting %s users (%s new)" % (len(store), len(newusers))
    transaction.commit()
    
    # update profiles for new users as a separate transaction
    for flag, user, row in newusers:
        prof = user.get_profile()
        prof.score = int(row['Reputation'])
        type = row['UserTypeId']
        if type == '4':
            prof.type = const.USER_MODERATOR
        if type == '5':
            prof.type = const.USER_ADMIN
        prof.save()
        
    print "*** updating %s user profiles" % len(newusers)
    transaction.commit()
    
    # remap to users
    users = dict( (key, value[1]) for key, value in store.items() )
    return users

#@transaction.commit_manually
def insert_posts(fname, limit, users):
    "Inserts the posts"

    store = {}

    rows = xml_reader(fname, limit=limit)
    
    for (index, row) in enumerate(rows):
        userid  = row.get('OwnerUserId')
        if not userid: # user has been destroyed
            continue
        PostTypeId = row['PostTypeId']
        id      = row['Id']
        body    = row['Body']
        views   = row['ViewCount']
        creation_date = parse_time(row['CreationDate'])
        author  = users[userid] # this is the user field
        post, flag = models.Post.objects.get_or_create(author=author, views=views, creation_date=creation_date)
        post.save()
        store[id] = (flag, post)
    
    newposts = filter(lambda x: x[0], store.values() )
    print "*** inserting %s posts (%s new)" % (len(store), len(newposts))
    transaction.commit()
    
    # remap to contain only posts
    posts = dict( (key, value[1]) for key, value in store.items() )
    return posts

@transaction.commit_manually
def insert_post_revisions(fname, limit, users, posts):
    """
    Inserts post revisions. Also responsible for parsing out closed/deleted 
    states from the post history log
    """
        
    rows = xml_reader(fname)
    
    # Stack overflow splits up modifications to title, tags, and content
    # as separate rows in the post history XML distinguished by type
    # We need to first go through and collect them together
    # based on the GUID they set on them.
    revisions = {} # Dictionary for fast GUID lookup
    guid_list = [] # Ensures order doesn't get mixed up
    for (index, row) in enumerate(rows):
        try:
            post   = posts[ row['PostId'] ]
            author = users[ row['UserId'] ]
        except KeyError:
            # post or author missing 
            continue
        
        guid    = row['RevisionGUID']
        datestr = row['CreationDate']
        date    = parse_time(datestr)
        
        if guid not in revisions:
            guid_list.append(guid)
        
        # this is the new revision we need to create
        # in our model a single revision contains all changes that were applied
        rev = revisions.get(guid, {'post':post, 'author':author, 'date':date})
        
        type = row['PostHistoryTypeId']
        
        if type in ['1', '4']: # Title
            rev['title'] = row['Text']
        
        if type in ['2', '5']: # Body
            rev['content'] = row['Text']
        
        # adds the tag
        if type in ['3', '6']: # Tags
            rev['tag_string'] = parse_tag_string(row['Text'])
            
        # applies each action to every post
        if type in ['10','11','12','13']: # Moderator actions
            actions = {'10':const.REV_CLOSE, '11':const.REV_REOPEN,
                       '12':const.REV_DELETE, '13':const.REV_UNDELETE}
            # this is defined in the models
            post.moderator_action(actions[type], author, date)
            
        revisions[guid] = rev
    
    transaction.commit()
    
    # Now we can actually insert the revisions
    for guid in guid_list:
        data = revisions[guid]
        post = data['post']
        del data['post']
        post.create_revision(**data)
        
    transaction.commit()
    print "*** inserted %s post revisions" % len(guid_list) 

@transaction.commit_manually
def insert_questions(fname, limit, posts):
    "Inserts questions and answers"
    
    quest_map, answ_map = {}, {}
    rows = xml_reader(fname, limit=limit)
    
    # filter out posts that have been removed before
    rows = filter(lambda x: x['Id'] in posts, rows)

    # broken up into separate steps to allow manual transactions
    
    # subselect the questions
    questions = filter(lambda x: x['PostTypeId'] == '1', rows)
    
    for row in questions:
        Id = row['Id']
        post = posts[ Id ]
        
        # update question with the title of the post
        title = row["Title"]
        tags  = row['Tags']
        post.title = title
        post.save()
        
        quest, flag = models.Question.objects.get_or_create(post=post)
        quest_map[Id] = quest
        
        # add tags if it is a newly created question
        if flag: 
            tag_string = parse_tag_string(tags)
            post.set_tags(tag_string)
        
        
    print "*** inserting %s questions" % len(quest_map)
    transaction.commit()
    
    # subselect the answers
    answers = filter(lambda x: x['PostTypeId'] == '2', rows)
    for row in answers:
        Id   = row['Id']
        post = posts[ Id ] 
           
        quest = quest_map[row['ParentId']]
        answ, flag = models.Answer.objects.get_or_create(question=quest, post=post)
        # add title if it is a newly created answer
        if flag:
            # will set a title to be easier to find in admin
            post.title = "A: %s" % quest.post.title
            answ.save()
            post.save()
        answ_map[Id] = answ
    
    print "*** inserting %s answers" % len(answ_map)        
    transaction.commit()
    
@transaction.commit_manually
def insert_votes(fname, limit, users, posts):
    store = {}
    votes = xml_reader(fname)
    for row in votes:
        post = posts.get(row['PostId'])
        user = users.get(row['UserId'])
        addr = users.get(row['IPAddress'])
        VoteType = row['VoteTypeId']
        
        if post and user:
            if VoteType == '1':
                vote_type = const.VOTE_ACCEPT
            elif VoteType == '2':
                vote_type = const.VOTE_UP
            elif VoteType == '3':
                vote_type = const.VOTE_DOWN
            else:
                continue
            
            vote, flag = models.Vote.objects.get_or_create(post=post, author=user, type=vote_type)
            store[row['Id']] = vote

    transaction.commit()
    print "*** inserted %s votes" % len(store)
 
@transaction.commit_manually
def insert_comments(fname, posts, users, limit):
    comment_map = {}
    rows = xml_reader(fname, limit=limit)
    
    for (index, row) in enumerate(rows):
        Id = row['Id']
        try:
            parent = posts[ row['PostId'] ]
        except KeyError:
            continue # We haven't inserted this post
        text   = row['Text']
        userid = row['UserId']
        author = users[userid]
        creation_date = parse_time(row['CreationDate'])
        #post = comment_post_map[Id]
        post = models.Post.objects.create(author=author, creation_date=creation_date)
        post.content = 'Original body not yet available!'
        post.html = text
        post.save()
        comment, flag = models.Comment.objects.get_or_create(parent=parent, post=post)
        comment_map[Id] = comment
    transaction.commit()
    
    print "*** inserted %s comments" % len(comment_map)
    
@transaction.commit_manually
def insert_badges(fname, limit):
    "Inserts the badges"
    store = {}
    rows = xml_reader(fname, limit=limit)
    for (index, row) in enumerate(rows):
        Id = row['Id']
        type = row['Class']
        name   = row['Name']
        desc = row['Description']
        unique = row['Single'] == 'true'
        secret = row['Secret'] == 'true'
        type = {'3':models.BADGE_BRONZE, '2':models.BADGE_SILVER, '1':models.BADGE_GOLD}[type]
        badge = models.Badge.objects.create(name=name, type=type, description=desc, unique=unique, secret=secret)
        badge.save()
        store[Id] = badge
    transaction.commit()
    print "*** inserted %s badges" % len(store)
    return store

@transaction.commit_manually
def insert_awards(fname, users, badges, limit):
    "Inserts the badge awards"
    store = {}
    rows = xml_reader(fname)
    for (index, row) in enumerate(rows):
        Id = row['Id']
        try:
            user = users[row['UserId']]
            badge = badges[row['BadgeId']]
        except KeyError:
            continue
        date = parse_time(row['Date'])
        a = models.Award.objects.create(user=user, badge=badge, date=date)
        a.save()
        store[Id] = a
    transaction.commit()
    print "*** inserted %s badge awards" % len(store)
    

def execute(path, limit=None):
    """
    Executes the imports
    """
    
    # insert users into the database
    fname = join(path, 'Users.xml')
    users = insert_users(fname=fname, limit=limit)
    
    fname = join(path, 'Posts.xml')
    posts = insert_posts(fname=fname, limit=limit, users=users)
    
    fname = join(path, 'PostHistory.xml')
    revisions = insert_post_revisions(fname=fname, limit=limit, posts=posts, users=users)
    
    fname = join(path, 'Posts.xml')
    insert_questions(fname=fname, limit=limit, posts=posts)

    fname = join(path, 'Posts2Votes.xml')
    insert_votes(fname=fname, limit=limit, posts=posts, users=users)
    
    fname = join(path, 'PostComments.xml')
    insert_comments(fname=fname, posts=posts, users=users, limit=limit)
    
    fname = join(path, 'Badges.xml')
    badges = insert_badges(fname=fname, limit=limit)
    
    fname = join(path, 'Users2Badges.xml')
    insert_awards(fname=fname, users=users, badges=badges, limit=limit)
    
if __name__ =='__main__':
    import doctest, optparse
    
    # for debugging
    #sys.argv.extend( ["-p", "se0"] )
    
    # options for the program
    parser = optparse.OptionParser()
    parser.add_option("-p", "--path", dest="path", help="directory or zip archive containing a full biostar SE1 datadump")
    parser.add_option("-L", "--limit", dest="limit", help="limit to these many rows per file", default=None)
    (opts, args) = parser.parse_args()
    
    # stop execution if no parameters were specified
    if not opts.path:
        parser.print_help()
        sys.exit()
        
    # also run the doctests
    doctest.testmod()
    
    # call into the main program
    execute(path=opts.path, limit=opts.limit)
    
