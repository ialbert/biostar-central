
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


def copy_posts(post_uids={}):

    #roots = PostsPost.objects.filter(type__in=Post.TOP_LEVEL)
    source = PostsPost.objects.all()

    for post in source:
        new_post = Post.objects.get_all(uid=post.id).first()
        # Skip if the posts exists.
        if new_post:
            continue
        content = util.strip_tags(post.content)
        author = User.objects.filter(profile__uid=post.author_id).first()
        lastedit_user = User.objects.filter(profile__uid=post.lastedit_user_id).first()

        new_post = Post(uid=post.id, html=post.html, type=post.type,
                        subs_count=post.subs_count, lastedit_user=lastedit_user,
                        author=author, status=post.status, rank=post.rank,
                        lastedit_date=post.lastedit_date,
                        content=content, comment_count=post.comment_count,
                        title=post.title,
                        vote_count=post.vote_count, thread_votecount=post.thread_score,
                        creation_date=post.creation_date, tag_val=post.tag_val,
                        reply_count=post.reply_count, book_count=post.book_count,
                        view_count=post.view_count)
        # Store tree for later processing
        root_tree = post_uids.setdefault(post.root_id, {})
        root_tree = root_tree.setdefault(post.parent_id, []).append(post.id)
        post_uids.setdefault(post.root_id, root_tree)

        print(new_post.has_accepted, "IN GENERATOR")

        yield new_post

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
        # Skip if post or author do not exist
        if not (post and author):
            continue

        new_vote = Vote(post=post, author=author, type=vote.type, uid=util.get_uuid(5))
        if post.author != author:
            # Update the user reputation only if the author is different.
            Profile.objects.filter(user=post.author).update(score=F('score') + 1)

        # The thread vote count represents all votes in a thread
        Post.objects.get_all(pk=post.root_id).update(thread_votecount=F('thread_votecount') + 1)
        auth.trigger_vote(vote_type=vote.type, post=post, change=1)

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
        post_uids = {}

        #print(Post.objects.filter(has_accepted=True).all())
        #1/0
        Post.objects.bulk_create(objs=copy_posts(post_uids=post_uids), batch_size=20)

        for root_uid in post_uids:
            root = Post.objects.filter(uid=root_uid).first()
            thread = post_uids[root_uid]
            for parent_uid in thread:
                parent = Post.objects.filter(uid=parent_uid).first()
                for single in thread[parent_uid]:
                    Post.objects.filter(uid=single).update(parent=parent, root=root)
                    print(Post.objects.get_all(uid=single).first().has_accepted, "SINGLE")
                thread_query = Post.objects.filter(status=Post.OPEN, root=root)
                reply_count = thread_query.exclude(uid=parent.uid).filter(type=Post.ANSWER).count()
                thread_score = thread_query.exclude(uid=root.uid).count()
                Post.objects.filter(parent=parent).update(reply_count=reply_count)
                Post.objects.filter(root=root).update(thread_score=thread_score)

        Vote.objects.bulk_create(objs=copy_votes(), batch_size=20)

        return
