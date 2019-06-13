import logging
from django.db.models import F
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.shortcuts import reverse
from taggit.models import Tag
from .models import Post, Award, Subscription
from . import tasks, auth
from django.db.models import F, Q

logger = logging.getLogger("biostar")


# @receiver(post_save, sender=Award)
# def send_award_message(sender, instance, created, **kwargs):
#     """
#     Send message to users when they receive an award.
#     """
#     template = "messages/awards_created.md"
#     badge_url = reverse('badge_view', kwargs=dict(uid=instance.badge.uid))
#     context = dict(badge_url=badge_url, award=instance, post=instance.post)
#
#     if created:
#         # Send local message synchronously
#         tasks.create_messages(template=template, extra_context=context, rec_list=[instance.user])
#     return
#

@receiver(post_save, sender=Post)
def finalize_post(sender, instance, created, **kwargs):

    # Determine the root of the post.
    root = instance.root if instance.root is not None else instance

    # Make last editor of the post the first in the list of contributors.
    root.thread_users.remove(instance.lastedit_user)
    root.thread_users.add(instance.lastedit_user)

    # Add tags
    instance.tags.clear()
    tags = [Tag.objects.get_or_create(name=name)[0] for name in instance.parse_tags()]
    instance.tags.add(*tags)

    # Update las contributor last editor of the root
    Post.objects.filter(uid=root.uid).update(lastedit_user=instance.lastedit_user,
                                             last_contributor=instance.last_contributor)

    if created:
        # Make the Uid user friendly
        instance.uid = instance.uid or f"p{instance.pk}"

        # Set the titles
        if instance.parent and not instance.title:
            instance.title = instance.parent.title

        # Only comments may be added to a parent that is answer or comment.
        if instance.parent and instance.parent.type in (Post.ANSWER, Post.COMMENT):
            instance.type = Post.COMMENT

        # Set post type if it was left empty.
        if instance.type is None:
            instance.type = Post.COMMENT if instance.parent else Post.FORUM

        # This runs only once upon object creation.
        instance.title = instance.parent.title if instance.parent else instance.title

        # Default tags
        instance.tag_val = instance.tag_val or "tag1,tag2"

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
            instance.title = "%s: %s" % (instance.get_type_display()[0], instance.root.title[:80])
            
        # Save the instance.
        instance.save()

        # Update all fields that need to bypass the instance save.

        # Update the root reply count for non toplevel posts.
        if not instance.is_toplevel:
            Post.objects.filter(id=instance.root.id).update(reply_count=F("reply_count") + 1)

        # Update the root answer count
        if instance.type == Post.ANSWER:
            Post.objects.filter(id=instance.root.id).update(answer_count=F("answer_count") + 1)

        # Update root comment count
        if instance.type == Post.COMMENT:
            Post.objects.filter(id=instance.root.id).update(comment_count=F("comment_count") + 1)
            Post.objects.filter(pk=instance.parent.pk, is_toplevel=False).update(comment_count=F("comment_count") + 1)

        # Update the parent reply counts.
        if instance.parent != instance.root:
            Post.objects.filter(pk=instance.parent.pk).update(reply_count=F("reply_count") + 1)

        # Create subscription for the author to the root.
        auth.create_subscription(post=instance.root, user=instance.author)

        # Send subscription messages to all subscribed users.
        tasks.notify_followers.spool(post=instance, author=instance.author)
