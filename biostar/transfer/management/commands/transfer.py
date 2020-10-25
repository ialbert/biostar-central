import logging
import time
import re
import os
from itertools import count, islice
import html2text


from django.core.management.base import BaseCommand
from django.conf import settings

from taggit.models import Tag

from biostar.accounts.models import User, Profile
from biostar.forum import util, markdown
from biostar.forum.models import Post, Vote, Subscription, Badge, Award
from biostar.transfer.models import UsersUser, PostsPost, PostsVote, PostsSubscription, BadgesAward, UsersProfile

logger = logging.getLogger("engine")

LIMIT = None

def timer_func():
    """
    Prints progress on inserting elements.
    """

    last = time.time()

    def elapsed(msg):
        nonlocal last
        now = time.time()
        sec = round(now - last, 1)
        last = now
        print(f"{msg} in {sec} seconds")

    def progress(index, step=5000, msg=""):
        nonlocal last
        if index % step == 0:
            elapsed(f"... {index} {msg}")

    return elapsed, progress


def uid_from_context(context):
    """
    Parse post.uid from a link
    """
    # id pattern found in link
    pattern = r"[0-9]+"
    uid = re.search(pattern, context)
    uid = uid.group(0) if uid else ""
    return uid


def bulk_copy_users(limit):

    current = dict()
    old_users = UsersUser.objects.order_by("id")

    def gen_users():

        logger.info(f"Transferring users")

        elapsed, progress = timer_func()

        # Allow limiting the input
        stream = islice(zip(count(1), old_users), limit)
        for index, user in stream:
            progress(index=index,  msg="users")

            username = f"{user.name.replace(' ', '-')}-{user.id}"
            # Create user
            new_user = User(username=username, email=user.email, password=user.password,
                            is_active=user.is_active, is_superuser=user.is_admin, is_staff=user.is_staff)

            current[user.email] = new_user

            yield new_user

    def gen_profile():
        logger.info(f"Transferring profiles")

        elapsed, progress = timer_func()

        # Allow limiting the input
        stream = islice(zip(count(1), old_users), limit)
        for index, user in stream:
            progress(index, msg="profiles")
            text = util.strip_tags(user.profile.info)

            # The incoming users have weekly digest prefs as a default.
            profile = Profile(uid=user.id, user=current.get(user.email), name=user.name,
                              message_prefs=user.profile.message_prefs,
                              role=user.type, last_login=user.last_login, html=user.profile.info,
                              date_joined=user.profile.date_joined, location=user.profile.location,
                              website=user.profile.website, scholar=user.profile.scholar, text=text,
                              watched_tags=user.profile.watched_tags,score=user.score,
                              twitter=user.profile.twitter_id, my_tags=user.profile.my_tags,
                              digest_prefs=user.profile.digest_prefs, new_messages=user.new_messages)

            yield profile

    # Bulk create the users, then profile.
    elapsed, progress = timer_func()

    User.objects.bulk_create(objs=gen_users(), batch_size=10000)
    ucount = User.objects.all().count()
    elapsed(f"transferred {ucount} users")

    Profile.objects.bulk_create(objs=gen_profile(), batch_size=10000)
    pcount = Profile.objects.all().count()
    elapsed(f"transferred {pcount} profiles")


def bulk_copy_votes(limit):
    posts = PostsVote.objects.all()
    posts_map = {post.uid: post for post in Post.objects.all()}
    users_map = {user.email: user for user in User.objects.all()}

    def gen_votes():
        logger.info("Transferring Votes")

        stream = zip(count(1), posts)
        stream = islice(stream, limit)

        elapsed, progress = timer_func()
        for index, vote in stream:
            progress(index, msg="votes")
            post = posts_map.get(str(vote.post_id))
            author = users_map.get(vote.author.email)
            # Skip existing votes and incomplete post/author information.
            if not (post and author) or not post.root:
                continue

            vote = Vote(post=post, author=author, type=vote.type, uid=vote.id,
                        date=vote.date)
            yield vote
    # Bulk create the users, then profile.
    elapsed, progress = timer_func()
    Vote.objects.bulk_create(objs=gen_votes(), batch_size=1000)
    vcount = Vote.objects.all().count()
    elapsed(f"transferred {vcount} votes")


def decode(s):
    """
    Replace null bytes in a string
    """
    return s.replace('\x00', '').replace('\0', '').replace('\000', '')


def add_tags(delete=False):
    logger.info("Transferring tags")

    if delete:
        # Delete tags before going forward.
        Tag.objects.all().delete()

    for post in Post.objects.iterator():
        tags = [t.strip() for t in post.parse_tags()]
        tags = [Tag.objects.get_or_create(name=name)[0] for name in tags]
        post.tags.remove(*tags)
        post.tags.add(*tags)


def bulk_copy_posts(limit):
    relations = {}
    all_users = User.objects.order_by("id")
    users_set = {user.profile.uid: user for user in all_users}

    def gen_posts():
        logger.info("Transferring posts")

        posts = PostsPost.objects.order_by("id")
        elapsed, progress = timer_func()
        stream = zip(count(1), posts)
        stream = islice(stream, limit)

        for index, post in stream:
            progress(index, msg="posts")

            author = users_set.get(str(post.author_id))
            lastedit_user = users_set.get(str(post.lastedit_user_id))
            # Incomplete author information loaded or existing posts.
            if not (author and lastedit_user):
                continue

            is_toplevel = post.type in Post.TOP_LEVEL

            rank = post.lastedit_date.timestamp()
            # Convert the content to markdown if its html

            force_text = post.content.strip().startswith("<")
            if force_text:
                try:
                    content = html2text.html2text(post.content, bodywidth=0)
                except Exception as exc:
                    content = post.content
                    logger.error(f"Failed parsing post={post.id}.")

                html = markdown.parse(content, clean=False, escape=False)
            else:
                content = post.content
                html = post.html

            new_post = Post(uid=post.id, html=decode(html), type=post.type, is_toplevel=is_toplevel,
                            lastedit_user=lastedit_user, thread_votecount=post.thread_score,
                            author=author, status=post.status, rank=rank, accept_count=int(post.has_accepted),
                            lastedit_date=post.lastedit_date, book_count=post.book_count,
                            content=decode(content), title=decode(post.title), vote_count=post.vote_count,
                            creation_date=post.creation_date, tag_val=decode(post.tag_val),
                            view_count=post.view_count)

            # Store parent and root for every post.
            relations[str(new_post.uid)] = [str(post.root_id), str(post.parent_id)]
            yield new_post

    def set_counts():
        logger.info("Setting post counts")
        descendants = lambda p: Post.objects.filter(root=p) if p.is_toplevel else Post.objects.filter(parent=p)
        posts = {post: (descendants(post).exclude(id=post.id).filter(type=Post.ANSWER).count(),
                        descendants(post).exclude(id=post.id).filter(type=Post.COMMENT).count())
                 for post in Post.objects.all()}

        for post in posts:
            answer_count, comment_count = posts[post][0], posts[post][1]
            post.reply_count = answer_count + comment_count
            post.answer_count = answer_count
            post.comment_count = comment_count
            yield post

    def gen_updates():
        logger.info("Updating post relations")
        posts = {post.uid: post for post in Post.objects.all()}
        for pid in relations:
            root_uid, parent_uid = relations[pid][0], relations[pid][1]
            post = posts[pid]
            root = posts.get(root_uid)
            parent = posts.get(parent_uid)
            if not (root and parent):
                continue
            post.root = root
            post.parent = parent
            yield post

    def update_threadusers():
        logger.info("Updating thread users.")
        roots = Post.objects.filter(type__in=Post.TOP_LEVEL)
        # Get all authors belonging to descendants of root
        get_authors = lambda root: Post.objects.filter(root=root).values_list("author")
        # Create dict key by root and its contributors
        roots = {root: User.objects.filter(pk__in=get_authors(root=root)) for root in roots}
        for post in roots:
            users = roots[post]
            if post.root.thread_users.filter(pk__in=users):
                continue

            post.thread_users.add(*users)

    def gen_awards():
        logger.info("Transferring awards.")
        # Query the badges
        badges = {badge.name: badge for badge in Badge.objects.all()}
        awards = BadgesAward.objects.all()

        stream = zip(count(1), awards)
        stream = islice(stream, limit)
        elapsed, progress = timer_func()
        for index, award in stream:
            progress(index, msg="awards")
            badge = badges[award.badge.name]
            user = users_set.get(str(award.user_id))
            # Get post uid from context
            post_uid = uid_from_context(award.context)
            post = Post.objects.filter(uid=post_uid).first()
            # Bail when a user or post do not exist
            if not user or (post_uid and not post):
                continue
            award = Award(date=award.date, post=post, badge=badge, user=user, uid=award.id)

            yield award

    elapsed, progress = timer_func()
    Post.objects.bulk_create(objs=gen_posts(), batch_size=1000)
    pcount = Post.objects.all().count()
    elapsed(f"transferred {pcount} posts")

    Post.objects.bulk_update(objs=gen_updates(), fields=["root", "parent"],
                             batch_size=1000)
    update_threadusers()
    add_tags()

    elapsed(f"Updated {pcount} post threads and added tags.")

    Post.objects.bulk_update(objs=set_counts(), fields=["reply_count", "comment_count", "answer_count"],
                             batch_size=1000)
    elapsed(f"Set {pcount} post counts.")

    Award.objects.bulk_create(objs=gen_awards(), batch_size=10000)
    acount = Award.objects.all().count()
    elapsed(f"transferred {acount} awards")


def bulk_copy_subs(limit):

    users = {user.profile.uid: user for user in User.objects.all()}
    posts = {post.uid: post for post in Post.objects.all()}

    def generate():
        subs = PostsSubscription.objects.order_by('-date')
        logger.info("Copying subscriptions")
        elapsed, progress = timer_func()
        stream = zip(count(1), subs)
        stream = islice(stream, limit)

        for index, sub in stream:
            progress(index, msg="subscriptions")
            user = users.get(str(sub.user_id))
            post = posts.get(str(sub.post_id))

            # Skip incomplete data or subs already made
            if not (user and post):
                continue

            sub = Subscription(uid=sub.id, type=sub.type, user=user, post=post, date=sub.date)

            yield sub

    def update_counts():
        logger.info("Updating post subs_count")
        # Recompute subs_count for
        posts_map = ((post, len(post.subs.exclude(user=post.author))) for post in Post.objects.all())

        for post, value in posts_map:
            post.subs_count = value
            yield post

    elapsed, progress = timer_func()
    Subscription.objects.bulk_create(objs=generate(), batch_size=10000)
    scount = Subscription.objects.all().count()
    elapsed(f"transferred {scount} subscriptions")

    Post.objects.bulk_update(objs=update_counts(), fields=["subs_count"], batch_size=10000)
    elapsed(f"Updated {scount} subscription counts")

    #dcount = Digest.objects.all().count()
    #Digest.objects.bulk_create_posts(objs=gen_digests(), batch_size=10000)
    #elapsed(f"Updated {dcount} user digests")

    return


def update_scores():

    users = UsersUser.objects.all()
    elapsed, progress = timer_func()

    stream = zip(count(1), users)
    stream = islice(stream, None)
    for idx, u in stream:
        # Get current user
        progress(idx, msg="users", step=100)
        current = User.objects.filter(profile__uid=u.id).first()

        if current and current.profile.score < u.score:
            Profile.objects.filter(pk=current.profile.pk).update(score=u.score)
    return


def update_votes():
    votes = PostsVote.objects.all()
    elapsed, progress = timer_func()
    stream = zip(count(1), votes)
    stream = islice(stream, None)
    for idx, v in stream:
        progress(idx, msg="votes", step=100)
        current = Vote.objects.filter(uid=f"{v.id}").first()
        if current:
            Vote.objects.filter(pk=current.pk).update(date=v.date)
    return


def test():

    user = UsersUser.objects.filter(id=1542).first()
    user2 = User.objects.filter(profile__uid='1542').first()
    print(user.name, user.score * 10)
    print(user2.profile.score)
    #1/0
    post_ids = [123258, 123260]
    seen = []
    for p in post_ids:
        subs = PostsSubscription.objects.filter(post_id=p)
        print(len(subs), p)
        for sub in subs:
            # print(sub.post.author.profile.digest_prefs)
            # print(sub.id,
            #       sub.user.id,
            #       sub.user.email,
            #       sub.user.profile.digest_prefs,
            #       sub.post.author.email,
            #       sub.post.title,
            #       sub.user.profile.message_prefs, sub.type)

            if sub.user.id in seen:
                print(f"Already subbed to post={p}, {sub.user.email}")
                continue

            seen.append(sub.user.id)

            print('-' * 100)
    print(seen)
    print(len(seen))
    return


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        parser.add_argument('--posts', action="store_true", help="Transfer posts from source database to target.")
        parser.add_argument('--users', action="store_true", help="Transfer users from source database to target.")
        parser.add_argument('--votes', action="store_true", help="Transfer votes from source database to target.")
        parser.add_argument('--subs', action="store_true", help="Transfer subs from source database to target.")
        parser.add_argument('--limit', '-n', type=int, help="Transfer subs from source database to target.")
        parser.add_argument('--tags', action="store_true", help="Add the tags to database ")
        parser.add_argument('--scores', action="store_true", help="Update user scores. ")
        parser.add_argument('--update_votes', action="store_true", help="Update votes in the database ")

    def handle(self, *args, **options):

        load_posts = options["posts"]
        load_users = options["users"]
        load_votes = options["votes"]
        load_subs = options["subs"]
        load_tags = options['tags']
        limit = options.get("limit") or LIMIT
        scores = options['scores']
        update_v = options['update_votes']

        print(f"OLD_DATABASE (source): {settings.OLD_DATABASE}")
        print(f"NEW_DATABASE (target): {settings.NEW_DATABASE}")

        #test()
        #return

        if load_posts:
            bulk_copy_posts(limit=limit)
            return
        if load_votes:
            bulk_copy_votes(limit=limit)
            return
        if load_users:
            bulk_copy_users(limit=limit)
            return
        if load_subs:
            bulk_copy_subs(limit=limit)
            return

        if load_tags:
            add_tags(delete=True)
            return

        if update_v:
            update_votes()
            return

        if scores:
            update_scores()
            return


        # Copy everything
        bulk_copy_users(limit=limit)

        bulk_copy_posts(limit=limit)

        bulk_copy_votes(limit=limit)

        bulk_copy_subs(limit=limit)

        return
