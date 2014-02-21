import sys, time, os, logging
import mailbox, markdown
from django.conf import settings
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User
from django.utils.timezone import utc
from django.utils import timezone, encoding
from email.utils import parsedate
from datetime import datetime
import itertools
from django.core.management.base import BaseCommand, CommandError
from optparse import make_option


logger = logging.getLogger(__name__)

def path_join(*args):
    return os.path.abspath(os.path.join(*args))


class Command(BaseCommand):
    help = 'migrate data from Biostar 1.*'

    option_list = BaseCommand.option_list + (
        make_option("-f", '--file', dest='file', default=False, help='import file'),
    )

    def handle(self, *args, **options):
        fname = options['file']
        if fname:
            summary(fname)

def parse_email(text):
    #print text
    address, name = text.split("(")[:2]
    address = address.replace(" at ", "@")
    address = address.replace(" ", "")
    address = address.lower()
    test = name.upper()
    if "ISO" in test or "UTF" in test:
        name = address.split("@")[0]
    name = name.strip(")")
    return address, name

class Bunch(object):
    pass

def no_junk(line):

    # invalid starts
    for word in "> --- From:".split():
        if line.strip().startswith(word):
            return False

    # junk words
    for word in "scrubbed attachment.html wrote: Sent:".split():
        if word in line:
            return False

    return True


def get_body(m):
    body = m.get_payload()
    if type(body) != str:
        body = '\n'.join([get_body(part.get_payload()) for part in body])

    # modify body to include at least one
    lines = body.splitlines()
    lines = filter(no_junk, lines)
    body = "\n".join(lines)
    body = markdown.markdown(body)
    return body

def valid(m):
    return m['Subject']

def unpack(m):
    b = Bunch()
    b.sender = m['From']
    b.id = m['Message-ID']
    b.reply_to = m['In-Reply-To']
    b.subj = m['Subject']
    assert b.subj, m
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

def create_post(b, author,  root=None, parent=None):
    title = b.subj
    body = b.body
    if not parent:
        title = title.strip()
        title = ' '.join(title.splitlines())
        title = ' '.join(title.split())
        title = title.title()
        post = Post(title=title, type=Post.QUESTION, content=body, tag_val="galaxy", author=author)
    else:
        post_type = Post.ANSWER if parent.is_toplevel else Post.COMMENT
        post = Post(type=post_type, content=body, tag_val="galaxy", author=author,  parent=parent)
    post.save()

    post.creation_date = post.lastedit_date = timezone.make_aware(b.datetime, timezone=utc)
    post.add_tags("galaxy")
    print "creating %s: %s" % (post.get_type_display(), title)
    return post

def summary(filename):



    mbox = mailbox.mbox(filename)
    mbox = filter(valid, mbox)
    mbox = map(unpack, mbox)

    mbox = filter(lambda m:  m.sender, mbox)
    mbox = map(fill, mbox)

    users = User.objects.all()
    users = dict( [(u.email, u) for u in users ])

    # create users as necessary
    for i, m in enumerate(mbox):
        if m.email not in users:
            print ("creating user %s" % m.email)
            username = m.email.replace("@",'.')[:30]
            u = User.objects.create(email=m.email)
            u.save()
            u.profile.name = m.display_name
            u.profile.save()
            users[m.email] = u


    # Get rid of posts
    Post.objects.all().delete()

    # need to build a tree of questions, answers and comments

    tree, posts = {}, {}

    for b in mbox:
        author = users[b.email]

        if not b.reply_to:
            post = create_post(b=b, author=author,)
            posts[b.id] = post
        else:
            parent = posts.get(b.reply_to)
            if parent:
                root = parent.root
                post = create_post(b=b, author=author, parent=parent)
                posts[b.id] = post



if __name__ == '__main__':
    #summary("mail/import.txt")
    #summary("mail/full.txt")
    pass
