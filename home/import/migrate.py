"""
Parses SE biostar xml archive, and populates a database with its
content then saves a loadable data fixture from this data.
"""
import sys, os, random, re
from datetime import datetime
from itertools import *
from xml.etree import ElementTree

# fixup path so we can access the Django settings module
join = os.path.join
file_dir = os.path.dirname(__file__)
sys.path.append( join(file_dir, '..' ))
sys.path.append( join(file_dir, '..', '..'))

# overide the default settings
os.environ['DJANGO_SETTINGS_MODULE'] = 'biostar_settings'

# test that the paths work
from django.conf import settings
from biostar.server import models
from django.db import transaction

def xml_reader(fname, limit=None):
    "XML dumps all use similar format, no attributes are used"
    elems = ElementTree.parse(fname)
    elems = islice(elems.findall('row'), limit)
    
    rows = []
    for elem in elems:
        # transforms children nodes to a dictionary keyed by tag with the node text as value"
        pairs = map( lambda x: (x.tag, x.text), elem)
        data  = dict(pairs)
        rows.append( data )

    return rows

def html_to_bbcode(body):
    body = body.replace("<", "[")
    body = body.replace(">", "]")
    return body

def parse_time(timestr):
    try:
        return datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        return datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S')

@transaction.commit_manually
def insert_users(fname, limit):
    "Inserts the users"
    store = {}
    rows = xml_reader(fname, limit=limit)
    for (index, row) in enumerate(rows):
        userid   = row['Id'] 
        username = row.get('Email', userid) or userid
        username = '%s%s' % (username, index)
        name = row.get('DisplayName', 'User %s' % userid)
        first_name = name
        try:
            u, f = models.User.objects.get_or_create(username=username, first_name=first_name)
        except Exception, e:
            print 'Failed inserting row %s' % row
            print userid, name, username
            raise(e)
        store[ userid ] = u
    transaction.commit()
    print "*** Inserted %s users" % len(store)
    return store

@transaction.commit_manually
def insert_posts(fname, user_map, limit):
    "Inserts the posts"
    store = {}
    rows = xml_reader(fname, limit=limit)
    for (index, row) in enumerate(rows):
        PostTypeId = row['PostTypeId']
        Id = row['Id']
        body   = row['Body']
        userid = row['OwnerUserId']
        views = row['ViewCount']
        creation_date = parse_time(row['CreationDate'])
        author = user_map[userid]
        body = html_to_bbcode(body)
        p = models.Post.objects.create(author=author, views=views, creation_date=creation_date)
        p.set(body)
        store[Id] = p
    transaction.commit()
    print "*** Inserted %s posts" % len(store)
    return store

@transaction.commit_manually
def insert_questions(fname, post_map, limit):
    "Inserts questions and answers"
    quest_map, answ_map = {}, {}
    rows = xml_reader(fname, limit=limit)
    
    tag_finder = re.compile('[^a-z]([a-z]+)[^a-z]')
    
    # broken up into separate steps to allow manual transactions
    for (index, row) in enumerate(rows):
        PostTypeId = row['PostTypeId']
        Id = row['Id']
        post = post_map[ Id ] 
        if PostTypeId == '1':
            # question
            title = row["Title"]
            tags = row['Tags']
            quest, flag = models.Question.objects.get_or_create(title=title, post=post)
            quest_map[Id] = quest
            if flag: # Only if newly created, add tags
                for tag in tag_finder.findall(tags):
                    tag = tag.strip()
                    if tag:
                        quest.tags.add(tag)
    
    transaction.commit()

    for (index, row) in enumerate(rows):
        PostTypeId = row['PostTypeId']
        Id = row['Id']
        post = post_map[ Id ] 
        if PostTypeId == '2':
            # answer
            quest = quest_map[row['ParentId']]
            answ, flag = models.Answer.objects.get_or_create(question=quest, post=post)
            answ_map[Id] = answ
    transaction.commit()
    
    print "*** Inserted %s questions" % len(quest_map)
    print "*** Inserted %s answers" % len(answ_map)
    
@transaction.commit_manually
def insert_votes(fname, user_map, post_map, limit):
    store = {}
    votes = xml_reader(fname)
    for row in votes:
        post = post_map.get(row['PostId'])
        user = user_map.get(row['UserId'])
        addr = user_map.get(row['IPAddress'])
        VoteType = row['VoteTypeId']
        # upmod=2, downmod=3
        valid = ('2', '3')
        
        if post and user and VoteType in valid:
            if VoteType == '2':
                vote_type = models.VOTE_UP
            else:
                vote_type = models.VOTE_DOWN
            
            v, flag = models.Vote.objects.get_or_create(post=post, author=user, type=vote_type)
            store[row['Id']] = v

    transaction.commit()
    print "*** Inserted %s votes" % len(store)
 
    
    
@transaction.commit_manually
def insert_comments(fname, post_map, user_map, limit):
    comment_map = {}
    rows = xml_reader(fname, limit=limit)
    
    for (index, row) in enumerate(rows):
        Id = row['Id']
        try:
            parent = post_map[ row['PostId'] ]
        except KeyError:
            continue # We haven't inserted this post
        text   = row['Text']
        userid = row['UserId']
        author = user_map[userid]
        creation_date = parse_time(row['CreationDate'])
        #post = comment_post_map[Id]
        post = models.Post.objects.create(author=author, creation_date=creation_date)
        post.set(text)
        comment, flag = models.Comment.objects.get_or_create(parent=parent, post=post)
        comment_map[Id] = comment
    transaction.commit()
    
    print "*** Inserted %s comments" % len(comment_map)
    
    

def execute(path, limit=300):
    """
    Executes the imports
    """
    fname = join(path, 'AnonUsers.xml')
    user_map = insert_users(fname=fname, limit=limit)

    fname = join(path, 'Posts.xml')
    post_map = insert_posts(fname=fname, limit=limit, user_map=user_map)

    fname = join(path, 'Posts.xml')
    insert_questions(fname=fname, limit=limit, post_map=post_map)

    fname = join(path, 'Posts2Votes.xml')
    insert_votes(fname=fname, limit=limit, post_map=post_map, user_map=user_map)
    
    fname = join(path, 'PostComments.xml')
    insert_comments(fname=fname, post_map=post_map, user_map=user_map, limit=limit)

if __name__ =='__main__':

    # requires a directory name as input
    # also runs from an editor (lot easier to work with initially)
    if len(sys.argv) == 1:
        dirname = 'datadump'
    else:
        dirname = sys.argv[1]

    execute(dirname)
    
