"""
Exports Biostar content into a tab delimited format
"""
__author__ = 'ialbert'
import sys, os
from main.server import models
pj = os.path.join

WDIR = "~/tmp/export"

WDIR = os.path.expanduser(WDIR)

def export_users(N):
    workdir = pj(WDIR, "about_me")

    limit = N or None

    fields = "id username display_name type uuid score website status location date_joined last_visited".split()

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    print "\t".join(fields)

    out = sys.stdout.write
    for user in models.User.objects.all().select_related("profile")[:limit]:
        p = user.profile
        fp = file(pj(workdir, str(user.id)), "wt")
        about_me = p.about_me or ''
        website = p.website or ''
        location = p.location or ''
        fp.write(about_me.encode("utf", "replace"))
        fp.close()
        dispay_name = p.display_name.encode("utf", "replace")

        data = [user.id, user.email, p.get_type_display(), p.uuid, p.display_name,
                p.score, p.scholar, p.my_tags, website, p.get_status_display(), location,
                user.date_joined, p.last_visited,
            ]

        data = map(unicode, data)
        line = u"\t".join(data)
        print line.encode("utf", "replace")

def export_posts(N):
    workdir = pj(WDIR, "posts")

    limit = N or None

    fields = "id root_id post_type parent.id author_id title tag_val views creation_date lastedit_date".split()

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    print "\t".join(fields)

    for post in models.Post.objects.all()[:limit]:

        try:
            title = post.title.encode("ascii", "replace")
            content = post.content.encode("ascii", "replace")
            tag_val = post.tag_val.encode("ascii", "replace")

            fp = file(pj(workdir, str(post.id)), "wt")
            fp.write(content)
            fp.close()

            data = [
                post.id, post.root.id, post.get_type_display(), post.parent.id, post.author.id, title,
                tag_val, post.views,
                post.creation_date, post.lastedit_date, post.lastedit_user.id
            ]
            data = map(unicode, data)
            print "\t".join(data)
        except Exception, e:
            sys.stderr.write("%s at %s\n" % (e, post.id))

def export_votes(N):

    limit = N or None

    fields = "author_id post_id vote_type vote_date".split()

    print "\t".join(fields)

    for vote in models.Vote.objects.all()[:limit]:

        data = [
            vote.author.id, vote.post.id, vote.get_type_display(), vote.date
            ]

        data = map(str, data)
        print "\t".join(data)

if __name__ == '__main__':
    import optparse
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-u", "--users", dest="users", action="store_true", help="prints the users to the standard out", default=0)
    parser.add_option("-p", "--posts", dest="posts", action="store_true", help="prints the posts to the standard out", default=0)
    parser.add_option("-v", "--votes", dest="votes", action="store_true", help="prints the votes to the standard out", default=0)

    parser.add_option("-n", dest="N", type=int, help="limits to N users", default=0)

    (opts, args) = parser.parse_args()

    if opts.users:
        export_users(opts.N)
        sys.exit()

    if opts.posts:
        export_posts(opts.N)
        sys.exit()

    if opts.votes:
        export_votes(opts.N)
        sys.exit()