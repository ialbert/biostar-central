import sys, time
import mailbox
from django.conf import settings
from main.server import models, const
from email.utils import parsedate
from datetime import datetime

def parse_email(text):
    #print text
    address, name = text.split("(")[:2]
    address = address.replace(" at ", "@")
    address = address.replace(" ", "")
    test = name.upper()
    if "ISO" in test or "UTF" in test:
        name = address.split("@")[0]
    name = name.strip(")")
    return address, name

class Bunch(object):
    pass

def get_body(m):
    body = m.get_payload()
    if type(body) == str:
        return body
    else:
        return '\n'.join([get_body(part.get_payload()) for part in body])

def unpack(m):
    b = Bunch()
    b.sender = m['From']
    b.id = m['Message-ID']
    b.reply_to = m['In-Reply-To']
    b.subj = m['Subject']
    b.subj = b.subj.replace("[galaxy-user]","").strip()
    b.body = get_body(m)
    # Sun, 2 Dec 2012 18:01:21 +0000

    b.date = m['Date']
    #b.datetime = time.strptime(b.date,"%a, %d %b %Y %H:%M:%S +0000")
    b.datetime = parsedate(b.date)
    b.datetime = datetime(*b.datetime[:6])
    b.m = m
    return b

def fill(b):
    b.email, b.display_name = parse_email(b.sender)
    return b

def create_post(b, author, post_type, root=None, parent=None):
    title = b.subj
    body = b.body
    post = models.Post(title=title, type=post_type, content=body, tag_val="galaxy", author=author, root=root, parent=parent)
    post.save()
    post.creation_date = post.lastedit_date = b.datetime
    post.set_tags()
    print "creating %s: %s" % (post.get_type_display(), title)
    return post

def summary(filename):

    mbox = mailbox.mbox(filename)

    mbox = map(unpack, mbox)
    mbox = filter(lambda m:  m.sender, mbox)
    mbox = map(fill, mbox)

    users = models.User.objects.all()
    users = dict( [(u.email, u) for u in users ])


    # create users as necessary
    for i, m in enumerate(mbox):
        if m.email not in users:
            print "creating user %s" % m.email, m.display_name
            username = m.email.split("@")[0]
            u = models.User.objects.create(username=username, email=m.email)
            u.profile.display_name = m.display_name
            u.profile.save()
            users[m.email] = u

    # get rid of posts
    #models.Post.objects.all().delete()

    # need to build a tree of questions, answers and comments

    tree, posts = {}, {}

    for b in mbox:
        author = users[b.email]

        if not b.reply_to:
            post_type = const.POST_QUESTION
            post = create_post(b=b, author=author, post_type=post_type)
            posts[b.id] = post
        else:
            parent = posts.get(b.reply_to)
            if parent:
                if parent.type == const.POST_QUESTION:
                    post_type = const.POST_ANSWER
                else:
                    post_type = const.POST_COMMENT
                root = parent.root
                post = create_post(b=b, author=author, post_type=post_type, root=root, parent=parent)
                posts[b.id] = post

if __name__ == '__main__':
    summary("mail/import.txt")
    #summary("mail/full.txt")
