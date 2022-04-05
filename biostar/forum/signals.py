import logging
from django.db.models.signals import post_save, pre_save
from django.dispatch import receiver
from taggit.models import Tag
from django.db.models import F, Q
from biostar.accounts.models import Profile, Message, User
from biostar.forum.models import Post, Award, Subscription, SharedLink, Diff
from biostar.forum import tasks, auth, util


logger = logging.getLogger("engine")


@receiver(post_save, sender=Award)
def send_award_message(sender, instance, created, **kwargs):
    """
    Send message to users when they receive an award.
    """
    if created:
        template = "messages/awards_created.md"
        context = dict(award=instance)

        # Temporarily stop messages to high rep users.
        # if instance.user.profile.score < 1000:
        # Send local message
        tasks.create_messages(template=template, extra_context=context, user_ids=[instance.user.pk])

    return


@receiver(post_save, sender=SharedLink)
def send_herald_message(sender, instance, created, **kwargs):
    """
    Send message to users when they receive an award.
    """
    if created:
        template = "messages/shared_link.md"
        context = dict(shared=instance)

        # Let the user know we have received.
        tasks.create_messages(template=template, extra_context=context, user_ids=[instance.author.pk])
        logmsg = f"{instance.get_status_display().lower()} herald story {instance.url[:100]}"
        auth.db_logger(user=instance.author, text=logmsg)

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
        # Remove all 'lastedit user' flags for this user.
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

    # Label all posts by a spammer as 'spam'
    if instance.is_spammer:
        Post.objects.filter(author=instance.user).update(spam=Post.SPAM)


@receiver(post_save, sender=Post)
def finalize_post(sender, instance, created, **kwargs):

    # Determine the root of the post.
    root = instance.root if instance.root is not None else instance

    # Update last contributor, last editor, and last edit date to the thread
    Post.objects.filter(uid=root.uid).update(lastedit_user=instance.lastedit_user,
                                             lastedit_date=instance.lastedit_date)

    # Get newly created subscriptions since the last edit date.
    subs = Subscription.objects.filter(date__gte=instance.lastedit_date, post=instance.root)
    extra_context = dict()

    if created:
        # Make the uid user friendly
        instance.uid = instance.uid or f"9{instance.pk}"

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

        # Make the last editor first in the list of contributors
        # Done on post creation to avoid moderators being added for editing a post.
        instance.root.thread_users.remove(instance.lastedit_user)
        instance.root.thread_users.add(instance.lastedit_user)

        # Update this post rank on create and not every edit.
        instance.rank = instance.lastedit_date.timestamp()

        # Save the instance.
        instance.save()
        instance.update_parent_counts()

        # Bump the root rank when a new answer is added.
        if instance.is_answer:
            Post.objects.filter(uid=instance.root.uid).update(rank=util.now().timestamp())

        # Create subscription to the root.
        auth.create_subscription(post=instance.root, user=instance.author)

        # Get all subscribed users when a new post is created
        subs = Subscription.objects.filter(post=instance.root)

        # Notify users who are watching tags in this post
        tasks.notify_watched_tags.spool(uid=instance.uid, extra_context=extra_context)

        # Send out mailing list when post is created.
        tasks.mailing_list.spool(uid=instance.uid, extra_context=extra_context)

    # Set the tags on the instance.
    if instance.is_toplevel:
        tags = [Tag.objects.get_or_create(name=name)[0] for name in instance.parse_tags()]
        instance.tags.clear()
        instance.tags.add(*tags)

    # Ensure spam posts get closed status
    if instance.is_spam:
        Post.objects.filter(uid=instance.uid).update(status=Post.CLOSED)

    if not instance.is_toplevel:
        # Title is inherited from top level.
        title = f"{instance.get_type_display()}: {instance.root.title[:80]}"
        Post.objects.filter(uid=instance.uid).update(title=title)

    # Ensure posts get re-indexed after being edited.
    Post.objects.filter(uid=instance.uid).update(indexed=False)

    # Exclude current authors from receiving messages from themselves
    subs = subs.exclude(Q(type=Subscription.NO_MESSAGES) | Q(user=instance.author))

    sub_ids = list(subs.values_list('id', flat=True))

    # Notify subscribers
    tasks.notify_followers.spool(sub_ids=sub_ids,
                                 author_id=instance.author.pk,
                                 uid=instance.uid,
                                 extra_context=extra_context)


@receiver(post_save, sender=Post)
def check_spam(sender, instance, created, **kwargs):
    # Classify post as spam/ham.
    tasks.spam_check.spool(uid=instance.uid)

@receiver(post_save, sender=SharedLink)
def link_title(sender, instance, created, **kwargs):
    # Set the title of each link upon creation
    if created:
        tasks.set_link_title.spool(pk=instance.pk)