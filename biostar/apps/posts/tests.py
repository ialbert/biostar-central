"""
Post related tests.

These will execute when you run "manage.py test".
"""

import logging
from django.conf import settings
from biostar.apps.people.models import User, Profile
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.notifications.models import Message

from django.test import TestCase

logging.disable(logging.CRITICAL)


class PostTest(TestCase):

    def test_post_creation(self):
        "Testing post creation"
        eq = self.assertEqual

        # Create an admin user and a post.
        title = "Hello World!"
        email = "john@this.edu"
        jane = User.objects.create(email=email)
        post = Post(title=title, author=jane, type=Post.FORUM)
        post.save()
        eq(post.type, Post.FORUM)
        eq(post.root, post)
        eq(post.parent, post)

        # Subscriptions are automatically created
        sub = Subscription.objects.get(user=jane)
        eq(sub.user, jane)
        eq(sub.post, post)

        # A new post triggers a message to the author.
        email = "jane@this.edu"
        john = User.objects.create(email=email)
        answer = Post(author=john, parent=post)
        answer.save()

        eq(answer.root, post)
        eq(answer.parent, post)
        eq(answer.type, Post.ANSWER)

        # Add comment. The parent will override the post type.
        email = "bob@this.edu"
        bob = User.objects.create(email=email)
        comment = Post(author=bob, type=Post.FORUM, parent=answer)
        comment.save()

        eq(comment.root, post)
        eq(comment.parent, answer)
        eq(comment.type, Post.COMMENT)

        # Everyone posting in a thread gets a subscription to the root post of the
        subs = Subscription.objects.filter(post=post)
        eq(len(subs), 3)
