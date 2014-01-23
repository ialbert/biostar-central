"""
Notification related tests.

These will execute when you run "manage.py test".
"""

import logging
from django.conf import settings
from biostar.apps.people.models import User, Profile
from biostar.apps.posts.models import Post, Subscription
from biostar.apps.notes.models import Message

from django.test import TestCase

logging.disable(logging.CRITICAL)


class NoteTest(TestCase):

    def test_note_creation(self):
        "Testing notifications"

        eq = self.assertEqual

        # Create some users
        title = "Hello Notifications!"
        emails = ["john@this.edu", "jane@this.edu", "bob@this.edu", "alice@this.edu",
                  "bill@this.edu", "jeff@this.edu" ]
        users, posts = [], []
        parent = None
        for email in emails:
            user = User.objects.create(email=email)
            users.append(user)

            post = Post(title=title, author=user, type=Post.FORUM, parent=parent)
            post.save()
            posts.append(post)

            parent = post

        count = len(emails)

        allp = Post.objects.all()
        eq(len(allp), count)

        answ = Post.objects.filter(type=Post.ANSWER)
        eq(len(answ), 1)

        comm = Post.objects.filter(type=Post.COMMENT)
        eq(len(comm), count - 2)

        subs = Subscription.objects.get_subs(post=parent)
        eq(len(subs), count)

        # now test the number of messages that they have
        for index, email in enumerate(emails):
            num = Message.objects.filter(user__email=email).count()
            eq(num, count - index - 1)



