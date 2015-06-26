from __future__ import absolute_import, division, unicode_literals

import logging

from django.test import TestCase
from django.conf import settings
from django.core import mail

from biostar3.forum import apps, models, auth
from biostar3.forum.models import User, Post, PostSub, Message

from biostar3.forum import html

from faker import Factory
from django.contrib.sites.models import Site

logging.disable(logging.INFO)

class Bunch(object):
    def __init__(self):
        self.sub_count, self.msg_count, self.email_count = -1, -1, -1

def snapshot(user=None):
    b = Bunch()

    if user:
        sub_count = PostSub.objects.filter(user=user).count()
        msg_count = Message.objects.filter(user=user).count()
    else:
        sub_count = PostSub.objects.count()
        msg_count = Message.objects.count()

    b.email_count = len(mail.outbox)
    b.sub_count = sub_count
    b.msg_count = msg_count

    return b

class SimpleTests(TestCase):
    def test_admin_user(self):
        """
        Test that admin users are created.
        """
        TRUE = self.assertTrue

        TRUE(len(settings.ADMINS) > 0)
        for user, email in settings.ADMINS:
            user = User.objects.get(email=email)
            TRUE(user.is_admin)

    def test_user_gen(self):
        """
        Tests that user welcome emails are sent.
        """
        EQ = self.assertEqual

        f = Factory.create()
        count = 10

        before = len(mail.outbox)

        for i in range(count):
            user = User.objects.create(name=f.name(), email=f.email())

        after = len(mail.outbox)

        # Test sending automated emails.
        EQ(before + count, after)

    def test_local_links(self):
        "Links to local content are reformatted"
        EQ = self.assertEqual
        f = Factory.create()
        site = Site.objects.get_current()
        user = User.objects.create(email=f.email(), name=f.name())

        user_link = "http://%s/u/%s" % (site.domain, user.id)
        user_text = "ABC %s ABC" % user_link
        user_html = '<p>ABC <a href="%s">%s</a> ABC</p>' % (user_link, user.name)

        pairs = [
            (user_text, user_html),
            ("a > b and `a > b`", "<p>a &gt; b and <code>a &gt; b</code></p>")
        ]

        for text, expect in pairs:
            result = html.sanitize(text, user).strip()
            EQ(expect, result)

    def test_embed_links(self):
        "Links to gist and youtube are embedded properly"
        EQ = self.assertEqual
        pairs = [
            # input, expected
            ("abcd", "abcd"),
            ("https://gist.github.com/123", "%s" % html.get_embedded_gist(123)),
            ("https://www.youtube.com/watch?v=123", "%s" % html.get_embedded_youtube(123)),
        ]
        for text, expected in pairs:
            result = html.embed_links(text)
            EQ(result, expected)



    def test_user_content_creation(self):
        "Creating a post subscribes the user to the post."

        EQ, TRUE = self.assertEqual, self.assertTrue

        f = Factory.create()
        jane = User.objects.create(email=f.email(), name=f.name())
        jane.profile.message_prefs = settings.SMART_MODE
        jane.profile.save()

        data = dict(
            title = f.sentence(),
            content= f.text(),
            tags = "hello world, goodbye",
        )

        before = snapshot()
        # Jane creates a post.
        post = auth.create_toplevel_post(user=jane, data=data)
        after1 = snapshot()

        # User has an EMAIL subscription but in SMART mode no email sent for this post.
        TRUE(PostSub.objects.filter(user=jane, post=post, type=settings.EMAIL_TRACKER))
        EQ(before.sub_count + 1, after1.sub_count)
        EQ(before.email_count, after1.email_count)
        EQ(before.msg_count, after1.msg_count)

        # Users on EMAIL_TRACKER will get emails.
        jane.profile.message_prefs = settings.EMAIL_TRACKER
        jane.profile.save()

        # User creates another toplevel post.
        post = auth.create_toplevel_post(user=jane, data=data)
        after2 = snapshot()

        # User has an EMAIL subscription but and an email was sent.
        TRUE(PostSub.objects.filter(user=jane, post=post, type=settings.EMAIL_TRACKER))
        EQ(before.sub_count + 2, after2.sub_count)
        EQ(before.email_count + 1, after2.email_count)
        EQ(before.msg_count, after2.msg_count)

        # A second user posts a followup to last post.
        joe = User.objects.create(email=f.email(), name=f.name())
        reply = auth.create_content_post(user=joe, parent=post, content=f.sentence())
        after3 = snapshot()

        # A new subscription but no email sent.
        EQ(before.sub_count + 3, after3.sub_count)
        EQ(before.email_count + 3, after3.email_count) # Also created a welcome email
        EQ(before.msg_count + 1, after3.msg_count)
