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
import textwrap

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

#
# Detecting the various way the From field may be formatted.
#
# User name followed by email in angled brackets: John Doe <foo@bar.com>
patt1 = OneOrMore(Word(chars)).setResultsName('name') + "<" + Word(email).setResultsName('email') + ">"

# An email alone: foo@bar.com
patt2 = Word(email).setResultsName('email')

# An email with angled brackerts: <foo@bar.com>
patt3 = '<' + Word(email).setResultsName('email') + '>'

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
    name = ' '.join(res.name)
    email = res.email

    # Corrects the name.
    name = check_name(name, email)
    return name, email


class Bunch(object):
    pass


def no_junk(line):
    "Gets rid of lines that contain junk"
    # invalid starts
    for word in "> --- From:".split():
        if line.strip().startswith(word):
            return False
    # junk words
    for word in "scrubbed attachment.html wrote: Sent:".split():
        if word in line:
            return False
    return True


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

    post.creation_date = post.lastedit_date = b.datetime
    post.save()

    #post.add_tags("galaxy")
    logger.info("creating %s: %s" % (post.get_type_display(), title))

    return post


REPLACE_PATT = [
    "[Galaxy-user]",
    "[Galaxy-User]",
    "[galaxy-user]",
    "[Bioc]",
]


def format_text(text):
    assert type(text), str
    lines = text.splitlines()
    lines = filter(no_junk, lines)
    lines = [textwrap.fill(line, width=100) for line in lines]
    text = "\n".join(lines)
    text = unicode(text, encoding="utf8", errors="replace")
    text = "<div class='preformatted'>" + text + "</div>"
    return text


def collect(m, data=[]):
    if m.is_multipart():
        for part in m.get_payload():
            data = collect(part, data)
    else:
        value = m.get_payload(decode=True)
        data.append(value)


def unpack_data(m):
    b = Bunch()
    b.name, b.email = parse_email(m)

    b.sender = m['From']
    b.id = m['Message-ID']
    b.reply_to = m['In-Reply-To']
    b.subj = m['Subject']

    b.subj = unicode(b.subj, encoding="utf8", errors="replace")
    assert b.subj, m

    for patt in REPLACE_PATT:
        b.subj = b.subj.replace(patt, '')

    try:
        data = []
        collect(m, data)
        b.body = "\n".join(data)
        b.body = format_text(b.body)

    except KeyError, exc:
        print exc
        print "skipping post %s" % b.subj


    b.date = m['Date']
    #b.datetime = time.strptime(b.date,"%a, %d %b %Y %H:%M:%S +0000")
    b.datetime = parsedate(b.date)
    b.datetime = datetime(*b.datetime[:6])

    b.datetime = timezone.make_aware(b.datetime, timezone=utc)

    return b

from django.db.models import signals

def parse_mbox(filename):

    signals.post_save.disconnect(dispatch_uid="create_messages")

    logger.info ("*** parsing mbox %s" % filename)

    # Parse the mbox.
    rows = mailbox.mbox(filename)

    # Keep only messages with a valid subject.
    mbox = ifilter(lambda m: m['Subject'], rows)

    # Create users
    rows = imap(unpack_data, rows)

    users = User.objects.all()
    users = dict([(u.email, u) for u in users])

    Post.objects.all().delete()

    logger.info("*** found %s users" % len(users))

    tree, posts = {}, {}

    #rows = islice(rows, 100)

    for b in rows:

        if b.email not in users:
            logger.info("--- creating user %s, %s" % (b.name, b.email))
            u = User.objects.create(email=b.email, name=b.name)
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

    logger.info("*** %s users" % len(users))
    logger.info("*** %s posts" % len(posts))


if __name__ == '__main__':
    pass
