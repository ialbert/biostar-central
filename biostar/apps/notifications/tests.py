"""
Notification related tests.

These will execute when you run "manage.py test".
"""

import logging
from django.conf import settings
from biostar.apps.people.models import User, Profile
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.notifications.models import Message

from django.test import TestCase

logging.disable(logging.CRITICAL)


class NoteTest(TestCase):

    def test_note_creation(self):
        "Testing notifications"

        eq = self.assertEqual

        # Create some users
        title = "Hello World!"
        emails = ["john@this.edu", "jane@this.edu", "bob@this.edu", "alice@this.edu", "bill@this.edu"]
        parent = None
        for email in emails:
            user = User.objects.create(email=email)
            post = Post(title=title, author=user, type=Post.FORUM, parent=parent)
            post.save()
            parent = post

        posts = Post.objects.all()
        eq(len(posts), 5)

        answ = Post.objects.filter(type=Post.ANSWER)
        eq(len(answ), 1)

        comm = Post.objects.filter(type=Post.COMMENT)
        eq(len(comm), 3)

        subs = Subscription.objects.get_subs(post=parent)
        eq(len(subs), 5)


