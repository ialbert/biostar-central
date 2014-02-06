"""
Exports Biostar content into a tab delimited format
"""
__author__ = 'ialbert'
import sys, os, shutil
from main.server import models
pj = os.path.join

MIGRATE_DIR = os.environ.get('BIOSTAR_MIGRATE_DIR')

if not MIGRATE_DIR:
    raise Exception("set the BIOSTAR_MIGRATE_DIR environment variable")

MIGRATE_DIR = os.path.expanduser(MIGRATE_DIR)

if not os.path.isdir(MIGRATE_DIR):
    os.makedirs(MIGRATE_DIR)

def export_users(N):
    workdir = pj(MIGRATE_DIR, "about_me")
    out_name = pj(MIGRATE_DIR, "users.txt")
    out_stream = file(out_name, 'wt')
    def write(line):
        out_stream.write('%s\n' % line)

    limit = N or None

    fields = "id type status email display_name score scholar my_tags website location date_joined last_visited".split()

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    write("\t".join(fields))

    for user in models.User.objects.all().select_related("profile")[:limit]:
        p = user.profile
        fp = file(pj(workdir, str(user.id)), "wt")
        about_me = p.about_me or ''
        website = p.website or ''
        location = p.location or ''
        fp.write(about_me.encode("utf", "replace"))
        fp.close()
        dispay_name = p.display_name.encode("utf", "replace")

        data = [user.id, p.get_type_display(), p.get_status_display(), user.email, p.display_name,
                p.score, p.scholar, p.my_tags, website,  location,
                user.date_joined, p.last_visited,
            ]

        data = map(unicode, data)
        line = u"\t".join(data)
        line = line.encode("utf", "replace")
        write(line)

    print ("*** wrote users into %s" % out_name)

def export_posts(N):
    workdir = pj(MIGRATE_DIR, "posts")

    out_name = pj(MIGRATE_DIR, "posts.txt")
    out_stream = file(out_name, 'wt')
    def write(line):
        out_stream.write('%s\n' % line)

    limit = N or None

    fields = "id root_id parent.id author_id post_type post_status title tag_val " \
        "score full_score views answer_count book_count accepted creation_date lastedit_date".split()

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    write("\t".join(fields))

    for post in models.Post.objects.all().order_by('id')[:limit]:

        try:
            title = post.title.encode("utf-8", "replace")
            html = post.html.encode("utf-8", "replace")
            tag_val = post.tag_val.encode("utf-8", "replace")

            fp = file(pj(workdir, str(post.id)), "wt")
            fp.write(html)
            fp.close()

            data = [
                post.id, post.root.id,  post.parent.id, post.author.id,
                post.get_type_display(), post.get_status_display(), title, tag_val,
                post.score, post.full_score, post.views, post.answer_count, post.book_count, post.accepted,
                post.creation_date, post.lastedit_date
            ]

            data = map(unicode, data)
            write("\t".join(data))

        except Exception, e:
            sys.stderr.write("%s at %s\n" % (e, post.id))

    print ("*** wrote posts into %s" % out_name)

def export_votes(N):

    out_name = pj(MIGRATE_DIR, "votes.txt")
    out_stream = file(out_name, 'wt')
    def write(line):
        out_stream.write('%s\n' % line)

    limit = N or None

    fields = "author_id post_id vote_type vote_date".split()

    write("\t".join(fields))

    for vote in models.Vote.objects.all()[:limit]:

        data = [
            vote.author.id, vote.post.id, vote.get_type_display(), vote.date
            ]

        data = map(str, data)
        write("\t".join(data))

    print ("*** wrote votes into %s" % out_name)

if __name__ == '__main__':
    import optparse
    usage = """usage: %prog [options]

BIOSTAR_MIGRATE_DIR=""" + MIGRATE_DIR

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--users", dest="users", action="store_true", help="prints the users to the standard out", default=0)
    parser.add_option("-p", "--posts", dest="posts", action="store_true", help="prints the posts to the standard out", default=0)
    parser.add_option("-v", "--votes", dest="votes", action="store_true", help="prints the votes to the standard out", default=0)

    parser.add_option("-n", dest="N", type=int, help="limits to N users", default=0)

    (opts, args) = parser.parse_args()

    print ("migration work directory: %s" % MIGRATE_DIR)

    if opts.users:
        export_users(opts.N)

    if opts.posts:
        export_posts(opts.N)

    if opts.votes:
        export_votes(opts.N)


