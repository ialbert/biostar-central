
import logging
from django.core.management.base import BaseCommand
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote
from biostar.forum import util
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
        logger.info(f"Created user email={new_user.email}")


def generate_posts():
    source = PostsPost.objects.all()
    for post in source:
        new_post = Post.objects.get_all(uid=post.id).first()
        # Skip if the posts exists.
        if new_post:
            continue
        content = util.strip_tags(post.content)

        yield Post(uid=post.id, html=post.html, type=post.type,
                   subs_count=post.subs_count, lastedit_user__profile_uid=post.lastedit_user_id,
                   author__profile_uid=post.author_id, root_id=post.root_id, status=post.status,
                   parent_id=post.parent_id, lastedit_date=post.lastedit_date,
                   content=content, comment_count=post.comment_count,
                   has_accepted=post.has_accepted, title=post.title,
                   thread_score=post.thread_score, vote_count=post.vote_count,
                   creation_date=post.creation_date, tag_val=post.tag_val,
                   reply_count=post.reply_count, book_count=post.book_count,
                   view_count=post.view_count)


def generate_votes():
    source = PostsVote.objects.all()
    for vote in source:
        vote


def copy_forum():
    """
    Bulk create posts from source database
    """
    # Copy users first
    copy_users()

    # Bulk create posts, votes, then subs in order.
    Post.objects.bulk_create(objs=generate_posts(), batch_size=20)
    Vote.objects.bulk_create(objs=generate_votes(), batch_size=20)
    return


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        # Get users from default database
        # Copy users, posts, votes, then subscriptions in order.
        copy_forum()


        return
