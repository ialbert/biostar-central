
import logging
from django.core.management.base import BaseCommand
from django.db.models import F
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote
from biostar.forum import util, auth
logger = logging.getLogger("engine")


def copy_users():

    source = UsersUser.objects.all()

    for user in source:

        new_user = User.objects.filter(email=user.email).first()
        if new_user:
            continue

        username = f"{user.name}{user.id}" if User.objects.filter(username=user.name).exists() else user.name

        # Create user
        new_user = User.objects.create(username=username, email=user.email,
                                       password=user.password, is_active=user.is_active,
                                       is_superuser=user.is_admin, is_staff=user.is_staff)
        # Update profile
        Profile.objects.filter(user=new_user).update(uid=user.id, name=user.name,
                              role=user.type, last_login=user.last_login,
                              date_joined=user.profile.date_joined, location=user.profile.location,
                              website=user.profile.website, scholar=user.profile.scholar, text=user.profile.info,
                              score=user.score, twitter=user.profile.twitter_id, my_tags=user.profile.my_tags,
                              digest_prefs=user.profile.digest_prefs, new_messages=user.new_messages)

    logger.info("Copied all users from biostar2")


def copy_posts(relations={}):
    logger.info("Copying posts.")
    source = PostsPost.objects.all().order_by("pk")[:10000]

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


# def update_posts(relations):
#
#     # Walk through forum tree and set the root, parent relations.
#     #
#     # def generate():
#     #
#     #     for post_uid in relations:
#     #         root_uid, parent_uid = relations[post_uid][0], relations[post_uid][1]
#     #         post = Post.objects.filter(uid=post_uid).first()
#     #         root = Post.objects.filter(uid=root_uid).first()
#     #         parent = Post.objects.filter(uid=parent_uid).first()
#     #
#     #         if not (root and parent):
#     #             continue
#     #         post.root = root
#     #         post.parent = parent
#     #         yield post
#     #
#     # Post.objects.bulk_update(objs=generate(), fields=["root", "parent"], batch_size=100)
#     return


def copy_votes():

    source = PostsVote.objects.all()[:1000]
    for vote in source:
        author = User.objects.filter(profile__uid=vote.author_id).first()
        post = Post.objects.get_all(uid=vote.post_id).first()
        new_vote = Vote.objects.filter(author=author, post=post, type=vote.type).first()
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

    for post_uid in relations:
        root_uid, parent_uid = relations[post_uid][0], relations[post_uid][1]
        root = Post.objects.filter(uid=root_uid).first()
        parent = Post.objects.filter(uid=parent_uid).first()
        Post.objects.filter(uid=post_uid).update(parent=parent, root=root)

        # Parent and root do not exist
        if not (root and parent):
            continue

        post = Post.objects.filter(uid=post_uid).first()
        thread = root if root is not None else post
        thread.thread_users.remove(post.author)
        thread.thread_users.add(post.author)
        thread_query = Post.objects.filter(status=Post.OPEN, root=root)

        reply_count = thread_query.exclude(uid=parent.uid).filter(type=Post.ANSWER).count()
        thread_score = thread_query.exclude(uid=root.uid).count()
        Post.objects.filter(parent=parent).update(reply_count=reply_count)
        Post.objects.filter(root=root).update(thread_score=thread_score)


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        #copy_users()
        copy_and_update()
        #Vote.objects.bulk_create(objs=copy_votes(), batch_size=20)

        return
