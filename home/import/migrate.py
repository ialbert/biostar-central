"""
Imports content from a directory that contains a SE biostar datadump
"""
import sys, os, random
from itertools import *
from xml.etree import ElementTree

join = os.path.join
curr_dir = os.path.dirname(__file__)
# fixup path so we can import directly
sys.path.append( join(curr_dir, '..', '..'))
sys.path.append( join(curr_dir, '..', '..', 'home'))

# overide the default settings
os.environ['DJANGO_SETTINGS_MODULE'] = 'biostar_settings'

import biostar
from django.conf import settings
from biostar.server import models

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

def execute(path, limit=100):
    """
    Imports data into the database
    """

    # users
    users = xml_reader(join(path, 'AnonUsers.xml'), limit=limit)
    user_map = {}
    for (index, row) in enumerate(users):
        #print row
        userid = row['Id'] 
        username = row.get('Email', userid) or userid
        username = '%s%s' % (username, index)
        name = row.get('DisplayName', 'User %s' % userid)
        first_name = name
        try:
            u, flag = models.User.objects.get_or_create(username=username, first_name=first_name)
        except Exception, e:
            print 'Failed inserting row %s' % row
            print userid, name, username
            raise(e)

        user_map[ userid ] = u
        
    print "*** Inserted %s users" % len(user_map)
    # migrate posts
    #
    posts = xml_reader(join(path, 'Posts.xml'), limit=limit)
    posts_map = {}
    
    for row in posts:
        post_type = row['PostTypeId']
        assert post_type in ('1', '2')
        # check for spa,
        if 'DeletionDate' in row:
            continue
        body = row['Body']
        userid = row['OwnerUserId']
        author = user_map[userid]


        # this will need to be transformed to BBcode but for now transform to HTML
        body = body.replace("&lt;", "<")
        body = body.replace("&gt;", ">")

        #print body 

        
        # create the post
        p, flag = models.Post.objects.get_or_create(html=body, author=author, lastedit_user=author)
        
        #print p, flag

        # create secondary entry
        if post_type == '1':
            # question
            title = row["Title"]
            q, flag = models.Question.objects.get_or_create(title=title, post=p)
            posts_map[row['Id']] = q

        elif post_type == '2':
            # answer
            q = posts_map[row['ParentId']]
            a, flag = models.Answer.objects.get_or_create(question=q, post=p)
            posts_map[row['Id']] = a
        else:
            # comment goes here
            pass

    print "*** Inserted %s posts" % len(posts_map)

    # add votes
    votes = xml_reader(join(path, 'Posts2Votes.xml'))
    user_rep = {}
    post_rep = {}

    for row in votes:
        post = posts_map.get(row['PostId'])
        user = user_map.get(row['UserId'])
        VoteType = row['VoteTypeId']
        # upmod=2, downmod=3
        valid = ('2', '3')
        
        if post and user and VoteType in valid:
            #print 'Creating %s' % row
            if VoteType == '2':
                vote = +10
            else:
                vote = -1
            user_rep.setdefault(user.id, []).append( vote )
            post_rep.setdefault(post.id, []).append( vote )
            v = models.Vote.objects.get_or_create(post=post.post, author=user, vote=vote)

    print "*** Inserted %s votes" % models.Vote.objects.all().count()
    
    # update the user and post reputations
    for userid, repdata in user_rep.items():
        user = models.User.objects.get(id=userid)
        user.score=sum(repdata)
        user.save()

    print "*** Updated %s user scores" % len(user_rep)
    
    for postid, repdata in post_rep.items():
        post = models.Post.objects.get(id=postid)
        post.views = random.randint(1, 100)
        post.votes = len(repdata)
        post.score = sum(repdata)
        post.save()
    
    print "*** Updated %s post scores" % len(user_rep)
    
    for question in models.Question.objects.all():
        question.answer_count = len(question.answers.all())
        question.save()
        
    print "*** Updated %s answer counts" % models.Question.objects.all().count()

if __name__ =='__main__':
    # requires a directory name as input
    dirname = sys.argv[1]
    #dirname = 'biostar-20110216151152'
    execute(dirname)