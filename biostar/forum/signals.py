import logging
from django.db.models.signals import post_save
from django.dispatch import receiver
from taggit.models import Tag
from django.db.models import F, Q
from biostar.accounts.models import Profile, Message, User
from .models import Post, Award, Subscription
from . import tasks, auth, util


logger = logging.getLogger("biostar")


@receiver(post_save, sender=Award)
def send_award_message(sender, instance, created, **kwargs):
    """
    Send message to users when they receive an award.
    """
    if created:
        template = "messages/awards_created.md"
        context = dict(award=instance)
        # Send local message
        tasks.create_messages(template=template, extra_context=context, rec_list=[instance.user])

    return


@receiver(post_save, sender=Profile)
def ban_user(sender, instance, created, **kwargs):
    """
    Delete all posts and awards belonging to a banned user.
    """

    if instance.state == Profile.BANNED:
        # Delete all posts by this users
        #print(Post.objects.filter(author=instance.user).thread_users)
        Post.objects.filter(author=instance.user).delete()
        #print(Post.objects.filter(author=instance))
        # Remove all 'lastedit user' flags
        # posts = Post.objects.filter(lastedit_user=instance.user)
        # for post in posts:
        #     Post.objects.filter(id=post.id).update(lastedit_user=post.author)

        # Delete all awards by the user.
        Award.objects.filter(user=instance.user).delete()

        Subscription.objects.filter(user=instance.user).delete()
        # Take out any personal information user added.
        #Profile.objects.filter(uid=instance.uid).update(text='')

        # Delete all messages
        Message.objects.filter(Q(sender=instance.user) | Q(recipient=instance.user)).delete()


@receiver(post_save, sender=Post)
def finalize_post(sender, instance, created, **kwargs):

    # Determine the root of the post.
    root = instance.root if instance.root is not None else instance

    # Add tags
    instance.tags.clear()
    tags = [Tag.objects.get_or_create(name=name)[0] for name in instance.parse_tags()]
    instance.tags.add(*tags)

    # Update last contributor, last editor, and last edit date to the thread
    Post.objects.filter(uid=root.uid).update(lastedit_user=instance.lastedit_user,
                                             last_contributor=instance.last_contributor,
                                             lastedit_date=instance.lastedit_date)

    # Get newly created subscriptions since the last edit date.
    subs = Subscription.objects.filter(date__gte=instance.lastedit_date, post=instance.root)

    if created:
        # Make the Uid user friendly
        instance.uid = instance.uid or f"p{instance.pk}"

        if instance.parent:
            # When the parent is set the root must follow the parent root.
            instance.root = instance.parent.root
        else:
            # When there is no parent, root and parent are set to itself.
            instance.root = instance.parent = instance

        # Answers and comments may only have comments associated with them.
        if instance.parent.type in (Post.ANSWER, Post.COMMENT):
            instance.type = Post.COMMENT

        # Sanity check.
        assert instance.root and instance.parent

        # Title is inherited from top level.
        if not instance.is_toplevel:
            instance.title = "%s: %s" % (instance.get_type_display(), instance.root.title[:80])

        # Make the last editor first in the list of contributors
        # Done on post creation to avoid moderators being added for editing a post.
        instance.root.thread_users.remove(instance.lastedit_user)
        instance.root.thread_users.add(instance.lastedit_user)

        # Update this post rank on create and not every edit.
        instance.rank = instance.lastedit_date.timestamp()

        # Save the instance.
        instance.save()
        instance.update_parent_counts()

        # Bump the root rank when a new descendant is added.
        Post.objects.filter(uid=instance.root.uid).update(rank=util.now().timestamp())

        # Create subscription to the root.
        auth.create_subscription(post=instance.root, user=instance.author)

        # Get all subscribed users when a new post is created
        subs = Subscription.objects.filter(post=instance.root)

        # Notify users who are watching tags in this post
        #tasks.notify_watched_tags(post=instance)

    # Ensure posts get re-indexed after being edited.
    Post.objects.filter(uid=instance.uid).update(indexed=False)

    # Exclude current authors from receiving messages from themselves
    subs = subs.exclude(Q(type=Subscription.NO_MESSAGES) | Q(user=instance.author))
    extra_context = dict(post=instance)
    tasks.notify_followers.spool(subs=subs, author=instance.author, extra_context=extra_context)
