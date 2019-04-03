
import logging
from django.core.management.base import BaseCommand
from django.db.models import F
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote
from biostar.forum import util, auth
logger = logging.getLogger("engine")


def copy_users():

    def create_users():
        source = UsersUser.objects.all()

        # Bulk create, then bulk update.
        for user in source:
            new_user = User.objects.filter(email=user.email).first()
            if new_user:
                continue
            username = f"{user.name}{user.id}"
            # Create user
            yield User(username=username, email=user.email, password=user.password,
                       is_active=user.is_active, is_superuser=user.is_admin, is_staff=user.is_staff)

    def create_profile():

        source = UsersUser.objects.all()
        for user in source:
            # Update profile
            profile = Profile.objects.filter(uid=user.id).first()
            if profile:
                continue
            new_user = User.objects.filter(email=user.email).first()
            yield Profile(uid=user.id, name=user.name, user=new_user,
                          role=user.type, last_login=user.last_login,
                          date_joined=user.profile.date_joined, location=user.profile.location,
                          website=user.profile.website, scholar=user.profile.scholar, text=user.profile.info,
                          score=user.score, twitter=user.profile.twitter_id, my_tags=user.profile.my_tags,
                          digest_prefs=user.profile.digest_prefs, new_messages=user.new_messages)

    # Bulk create the users, then profile.
    User.objects.bulk_create(objs=create_users(), batch_size=10000)
    Profile.objects.bulk_create(objs=create_profile(), batch_size=10000)

    logger.info("Copied all users from biostar2")


def copy_posts(relations={}):
    logger.info("Copying posts.")
    source = PostsPost.objects.all().order_by("pk")[:100000]

    for post in source:
        new_post = Post.objects.get_all(uid=post.id).first()
        # Skip if the posts exists.
        if new_post:
            continue
        content = util.strip_tags(post.content)
        author = User.objects.filter(profile__uid=post.author_id).first()
        lastedit_user = User.objects.filter(profile__uid=post.lastedit_user_id).first()

        new_post = Post(uid=post.id, html=post.html, type=post.type,
                        lastedit_user=lastedit_user,
                        author=author, status=post.status, rank=post.rank,
                        lastedit_date=post.lastedit_date,
                        content=content, title=post.title,
                        creation_date=post.creation_date, tag_val=post.tag_val,
                        view_count=post.view_count)

        # Store parent and root for every post.
        relations[new_post.uid] = [post.root_id, post.parent_id]
        yield new_post

    logger.info("Copied all posts from biostar2")


def copy_votes():

    source = PostsVote.objects.all()[:1000]
    for vote in source:
        author = User.objects.filter(profile__uid=vote.author_id).first()
        post = Post.objects.get_all(uid=vote.post_id).first()
        new_vote = Vote.objects.filter(uid=vote.id).first()
        # Skip if vote already exists.
        if new_vote:
            continue
        # Skip if post or author do not exist
        if not (post and author):
            continue

        new_vote = Vote(post=post, author=author, type=vote.type, uid=vote.id)
        if post.author != author:
            # Update the user reputation only if the author is different.
            Profile.objects.filter(user=post.author).update(score=F('score') + 1)

        # The thread vote count represents all votes in a thread
        auth.trigger_vote(vote_type=vote.type, post=post, change=+1)

        yield new_vote

    logger.info("Copied all votes from biostar2")


def copy_and_update():
    relations = {}
    Post.objects.bulk_create(objs=copy_posts(relations=relations), batch_size=1000)
    # Walk through tree and update parent, root, post, relationships
    logger.info("Updating posts")

    def generate_counts():

        for post_uid in relations:
            post = Post.objects.filter(uid=post_uid).first()
            post.root.thread_users.remove(post.author)
            post.root.thread_users.add(post.author)
            thread_query = Post.objects.filter(status=Post.OPEN, root=post.root)

            if post.parent == post:
                reply_count = thread_query.exclude(uid=post.parent.uid).filter(type=Post.ANSWER).count()
                post.reply_count = reply_count

            if post.is_toplevel:
                thread_score = thread_query.exclude(uid=post.uid).count()
                post.thread_score = thread_score

            yield post

    def generate():

        for post_uid in relations:
            root_uid, parent_uid = relations[post_uid][0], relations[post_uid][1]
            post = Post.objects.filter(uid=post_uid).first()
            root = Post.objects.filter(uid=root_uid).first()
            parent = Post.objects.filter(uid=parent_uid).first()

            if not root:
                continue
            post.root = root
            post.parent = parent
            yield post

    Post.objects.bulk_update(objs=generate(), fields=["root", "parent"], batch_size=100)
    Post.objects.bulk_update(objs=generate_counts(), fields=["thread_score", "reply_count"], batch_size=100)


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        copy_users()
        #

        #copy_and_update()
        #Vote.objects.bulk_create(objs=copy_votes(), batch_size=20)

        return
