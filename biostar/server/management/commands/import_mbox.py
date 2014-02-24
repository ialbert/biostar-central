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
from itertools import *

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
            parse_mbox(fname)


from pyparsing import Word, Or, alphanums, OneOrMore

chars = alphanums + r",?.-=_\/()[]"
email = alphanums + "+@._-"

patt1 = OneOrMore(Word( chars )).setResultsName('name') + "<" + Word( email ).setResultsName('email') + ">"
patt2 = Word( email ).setResultsName('email')
patt3 = '<' + Word( email ).setResultsName('email') + '>'

user_patt = patt1 | patt2 | patt3

def check_name(name, email):
    start = email.split("@")[0]
    if not name or '?' in name:
        return start

    return name.title()

def parse_email(data):
    sender = data['From']

    # print sender
    sender = sender.replace('"', '')
    check = sender.upper()
    res = user_patt.parseString(sender)
    name  = ' '.join(res.name)
    email = res.email

    # Corrects the name.
    name = check_name(name, email)
    return name, email



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

def format_text(text):
    assert type(text), str
    lines = text.splitlines()
    lines = filter(no_junk, lines)
    text = "\n".join(lines)
    text = unicode(text, encoding="utf8", errors="replace")
    text = markdown.markdown(text)
    return text

def get_body(m):

    if type(m) == str:
        body = format_text(m)
        return body

    if type(m) == list:
        body = '\n'.join([get_body(part.get_payload()) for part in m])
        return body

    body = m.get_payload()

    if type(body) == str:
        body = format_text(body)
        return body

    body = '\n'.join([get_body(part.get_payload()) for part in body])
    return body

def create_post(b, author, root=None, parent=None):
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
        post = Post(type=post_type, content=body, tag_val="galaxy", author=author, parent=parent)

    post.save()

    post.creation_date = post.lastedit_date = b.datetime
    post.add_tags("galaxy")
    print "creating %s: %s" % (post.get_type_display(), title)

    return post

REPLACE_PATT = [
    "[Galaxy-user]",
    "[Galaxy-User]",
    "[galaxy-user]"
]
def unpack_data(m):
    b = Bunch()
    b.name, b.email = parse_email(m)

    b.sender = m['From']
    b.id = m['Message-ID']
    b.reply_to = m['In-Reply-To']
    b.subj = m['Subject']
    assert b.subj, m

    for patt in REPLACE_PATT:
        b.subj = b.subj.replace(patt, '')

    try:
        b.body = get_body(m)
    except KeyError, exc:
        print exc
        print "skipping post %s" % b.subj
        print dir(m)
        sys.exit()
        b.body = ''

    b.date = m['Date']
    #b.datetime = time.strptime(b.date,"%a, %d %b %Y %H:%M:%S +0000")
    b.datetime = parsedate(b.date)
    b.datetime = datetime(*b.datetime[:6])

    b.datetime = timezone.make_aware(b.datetime, timezone=utc)

    return b


def parse_mbox(filename):
    print ("*** parsing mbox %s" % filename)

    # Parse the mbox.
    rows = mailbox.mbox(filename)

    # Keep only messages with a valid subject.
    mbox = ifilter(lambda m: m['Subject'], rows)

    # Create users
    rows = imap(unpack_data, rows)

    users = User.objects.all()
    users = dict([(u.email, u) for u in users])

    #Post.objects.all().delete()

    print "*** found %s users" % len(users)

    tree, posts = {}, {}

    rows = islice(rows, 10000)

    for b in rows:
        if b.email not in users:
            print ("--- creating user %s, %s" % (b.name, b.email))

            u = User.objects.create(email=b.email)
            u.save()
            u.profile.date_joined = b.datetime
            u.profile.save()
            users[u.email] = u

        #print b.body
        #print '-' * 20

        author = users[b.email]
        if not b.reply_to:
            post = create_post(b=b, author=author, )
            posts[b.id] = post
        else:
            parent = posts.get(b.reply_to)
            if parent:
                root = parent.root
                post = create_post(b=b, author=author, parent=parent)
                posts[b.id] = post


if __name__ == '__main__':

    pass
