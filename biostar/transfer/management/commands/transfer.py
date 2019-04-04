
import logging
from django.core.management.base import BaseCommand
from django.db.models import F
from biostar.accounts.models import User, Profile
from biostar.transfer.models import UsersUser, PostsPost, PostsVote
from biostar.forum.models import Post, Vote
from biostar.forum import util
logger = logging.getLogger("engine")


def bulk_copy_users():
    logger.info(f"Copying users")

    def gen_users():
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

        exists = list(User.objects.exclude(profile=None).values_list("email", flat=True))
        new_users = {user.email: user for user in User.objects.filter(profile=None)}
        source = UsersUser.objects.exclude(email__in=exists)

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


def copy_votes():

    source = PostsVote.objects.all()[:5000]

    def gen_votes():
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

            Post.objects.filter(uid=post.uid).update(vote_count=F('vote_count') + 1)

            if vote.type == Vote.BOOKMARK:
                Post.objects.filter(uid=post.uid).update(book_count=F('book_count') + 1)

            elif vote.type == Vote.ACCEPT:
                # There does not seem to be a negation operator for F objects.
                post.has_accepted = not post.has_accepted
                Post.objects.get_all(uid=post.root.uid).update(has_accepted=not post.root.has_accepted)

            Post.objects.get_all(root=post.root).update(thread_votecount=F('thread_votecount') + 1)
            Profile.objects.filter(user=post.author).update(score=F('score') + 1)

            yield Vote(post=post, author=author, type=vote.type, uid=vote.id)

    Vote.objects.bulk_create(objs=gen_votes(), batch_size=1000)

    logger.info("Copied all votes from biostar2")


def bulk_copy_and_update():
    relations = {}
    Post.objects.bulk_create(objs=copy_posts(relations=relations), batch_size=1000)
    # Walk through tree and update parent, root, post, relationships
    logger.info("Updating posts")

    def update_relationship():

        for post_uid in relations:
            root_uid, parent_uid = relations[post_uid][0], relations[post_uid][1]
            post = Post.objects.filter(uid=post_uid).first()
            root = Post.objects.filter(uid=root_uid).first()
            parent = Post.objects.filter(uid=parent_uid).first()

            if not root:
                continue
            root_thread = PostsPost.objects.filter(status=Post.OPEN, root_id=root_uid).values_list("id", flat=True).distinct()
            thread_query = Post.objects.filter(status=Post.OPEN, uid__in=list(root_thread))

            reply_count = thread_query.exclude(uid=post.uid).filter(type=Post.ANSWER).count()
            post.reply_count = reply_count

            post.thread_score = thread_query.exclude(uid=post.uid).count()

            post.root = root
            post.parent = parent
            post.root.thread_users.remove(post.author)
            post.root.thread_users.add(post.author)
            yield post

    Post.objects.bulk_update(objs=update_relationship(), fields=["root", "parent",  "thread_score", "reply_count"],
                             batch_size=1000)


class Command(BaseCommand):
    help = "Migrate users from one database to another."

    def add_arguments(self, parser):
        # Give the database file
        pass

    def handle(self, *args, **options):

        bulk_copy_users()
        #
        # copy_users()
        #bulk_copy_and_update()
        #copy_votes()

        return
