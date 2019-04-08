
import logging
from django.core.management.base import BaseCommand
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote
from biostar.forum import util
logger = logging.getLogger("engine")


def bulk_copy_users():

    def gen_users():
        logger.info(f"Copying users")
        exists = list(User.objects.values_list("email", flat=True))
        source = UsersUser.objects.exclude(email__in=exists)
        # Only iterate over users that do not exist already
        for user in source:
            username = f"{user.name}{user.id}"
            # Create user
            user = User(username=username, email=user.email, password=user.password,
                        is_active=user.is_active, is_superuser=user.is_admin, is_staff=user.is_staff)
            yield user

    def gen_profile():
        logger.info(f"Copying profiles")
        # Build user query
        new_users = {user.email: user for user in User.objects.filter(profile=None)}
        # Get users without profiles in target database.
        target = list(User.objects.exclude(profile=None).values_list("email", flat=True))
        # Exclude existing users from source database.
        source = UsersUser.objects.exclude(email__in=target)

        for user in source:
            profile = Profile(uid=user.id, name=user.name, user=new_users.get(user.email),
                              role=user.type, last_login=user.last_login,
                              date_joined=user.profile.date_joined, location=user.profile.location,
                              website=user.profile.website, scholar=user.profile.scholar, text=user.profile.info,
                              score=user.score, twitter=user.profile.twitter_id, my_tags=user.profile.my_tags,
                              digest_prefs=user.profile.digest_prefs, new_messages=user.new_messages)

            yield profile
    # Bulk create the users, then profile.
    User.objects.bulk_create(objs=gen_users(), batch_size=10000)
    Profile.objects.bulk_create(objs=gen_profile(), batch_size=10000)

    logger.info("Copied all users/profiles from biostar2")


def bulk_copy_votes():

    target = list(Vote.objects.values_list("uid", flat=True))
    source = PostsVote.objects.exclude(id__in=target)

    def gen_votes():
        posts = {post.uid: post for post in Post.objects.all()}
        users = {user.profile.uid: user for user in User.objects.all()}

        for vote in source:
            post = posts.get(vote.post_id)
            author = users.get(vote.author_id)
            # Skip if post or author do not exist
            if not (post and author):
                continue
            vote = Vote(post=post, author=author, type=vote.type, uid=vote.id)
            yield vote

    Vote.objects.bulk_create(objs=gen_votes(), batch_size=1000)

    logger.info("Copied all votes from biostar2")


def bulk_copy_posts():
    relations = {}

    # Walk through tree and update parent, root, post, relationships
    def gen_posts():
        logger.info("Copying posts.")
        exists = list(Post.objects.values_list("uid", flat=True))
        exists = [int(post_uid) for post_uid in exists if post_uid.isdigit()]

        source = PostsPost.objects.exclude(id__in=exists)
        users = {user.profile.uid: user for user in User.objects.all()}
        logger.info("Starting create loop.")
        for post in source:

            author = users[str(post.author_id)]
            lastedit_user = users[str(post.lastedit_user_id)]
            content = util.strip_tags(post.content)
            new_post = Post(uid=post.id, html=post.html, type=post.type,
                            lastedit_user=lastedit_user, thread_votecount=post.thread_score,
                            author=author, status=post.status, rank=post.rank, has_accepted=post.has_accepted,
                            lastedit_date=post.lastedit_date, book_count=post.book_count, reply_count=post.reply_count,
                            content=content, title=post.title,vote_count=post.vote_count, thread_score=post.reply_count,
                            creation_date=post.creation_date, tag_val=post.tag_val,
                            view_count=post.view_count)

            # Store parent and root for every post.
            relations[str(new_post.uid)] = [str(post.root_id), str(post.parent_id)]
            yield new_post

        logger.info("Copied all posts from biostar2")

    def gen_updates():
        logger.info("Updating post relations")
        posts = {post.uid: post for post in Post.objects.all()}
        for post_uid in relations:
            root_uid, parent_uid = relations[post_uid][0], relations[post_uid][1]
            post = posts[post_uid]
            root = posts.get(root_uid)
            parent = posts.get(parent_uid)
            if not (root and parent):
                continue
            post.root = root
            post.parent = parent
            yield post

        logger.info("Updated all posts from biostar2")

    def update_threadusers():
        logger.info("Updating thread users.")
        roots = Post.objects.filter(type__in=Post.TOP_LEVEL)
        roots = {root: Post.objects.filter(root=root) for root in roots}

        roots = {post: User.objects.filter(pk__in=roots[post].values_list("author")) for post in roots }

        for post in roots:
            users = roots[post]
            if post.root.thread_users.filter(pk__in=users):
                continue

            post.thread_users.add(*users)

    Post.objects.bulk_create(objs=gen_posts(), batch_size=1000)
    Post.objects.bulk_update(objs=gen_updates(), fields=["root", "parent"], batch_size=1000)
    update_threadusers()


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        #bulk_copy_users()
        #
        # copy_users()
        bulk_copy_posts()
        #bulk_copy_votes()

        return
