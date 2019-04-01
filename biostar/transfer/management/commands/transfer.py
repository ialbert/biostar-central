
import logging
from django.core.management.base import BaseCommand
from django.db.models import F
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote, trigger_vote
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


def copy_posts():

    source = PostsPost.objects.all().order_by("pk")

    for post in source:

        new_post = Post.objects.get_all(uid=post.id).first()
        # Skip if the posts exists.
        if new_post:

            continue
        content = util.strip_tags(post.content)
        author = User.objects.filter(profile__uid=post.author_id).first()
        lastedit_user = User.objects.filter(profile__uid=post.lastedit_user_id).first()
        parent = Post.objects.get_all(uid=post.parent_id).first()
        root = Post.objects.get_all(uid=post.root_id).first()

        Post.objects.create(uid=post.id, html=post.html, type=post.type,
                            subs_count=post.subs_count, lastedit_user=lastedit_user,
                            author=author, status=post.status, root=root, rank=post.rank,
                            parent=parent, lastedit_date=post.lastedit_date,
                            content=content, comment_count=post.comment_count,
                            has_accepted=post.has_accepted, title=post.title,
                            vote_count=post.vote_count, thread_votecount=post.thread_score,
                            creation_date=post.creation_date, tag_val=post.tag_val,
                            reply_count=post.reply_count, book_count=post.book_count,
                            view_count=post.view_count)
    logger.info("Copied all posts from biostar2")


def copy_votes():

    source = PostsVote.objects.all()
    for vote in source:
        author = User.objects.filter(profile__uid=vote.author_id).first()
        post = Post.objects.get_all(uid=vote.post_id).first()
        new_vote = Vote.objects.filter(author=author, post=post, type=vote.type).first()

        # Skip if vote already exists.
        if new_vote:
            continue

        new_vote = Vote(post=post, author=author, type=vote.type, uid=util.get_uuid(5))
        if post.author != author:
            # Update the user reputation only if the author is different.
            Profile.objects.filter(user=post.author).update(score=F('score') + 1)

        # The thread vote count represents all votes in a thread
        Post.objects.get_all(pk=post.root_id).update(thread_votecount=F('thread_votecount') + 1)
        trigger_vote(vote_type=vote.type, post=post, change=1)

        yield new_vote

    logger.info("Copied all votes from biostar2")


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        # Copy users first
        copy_users()
        copy_posts()
        Vote.objects.bulk_create(objs=copy_votes(), batch_size=20)

        return
