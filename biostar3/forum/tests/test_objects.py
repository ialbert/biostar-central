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


def create_user():
    f = Factory.create()
    user = User.objects.create(email=f.email(), name=f.name())
    user.profile.message_prefs = settings.SMART_MODE
    user.profile.save()
    return user


def create_post(user, parent=None, type=Post.QUESTION, title=None, tag_val=None, content=None):
    f = Factory.create()
    title = title or f.sentence()
    tag_val = tag_val or ", ".join(f.words())
    content = content or f.text()

    data = dict(
        title=title,
        content=content,
        tag_val=tag_val,
        type=type,
    )

    if not parent:
        post = auth.create_toplevel_post(user=user, data=data)
    else:
        post = auth.create_content_post(user=user, parent=parent, content=content)

    return post


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

    def setUp(self):
        self.jane = create_user()
        self.joe = create_user()
        self.EQ, self.TRUE = self.assertEqual, self.assertTrue

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
        EQ, TRUE = self.assertEqual, self.assertTrue

        # Bleach linkify will rewrite the html attributes.
        # So some outputs we cannot easily test for.

        # input, expected
        pairs = [
            ("https://gist.github.com/123", html.GIST_HTML.format(123)),

            # Bleach.linkfy rewrites attribute order.
            # Youtube embedding is tested only via a small fragment.for
            ('https://www.youtube.com/watch?v=1-2-3', 'src="//www.youtube.com/embed/1-2-3'),
            ('https://www.youtube.com/embed/X_Y-Z', 'src="//www.youtube.com/embed/X_Y-Z'),
            ('https://youtu.be/xyz', 'src="//www.youtube.com/embed/xyz'),

        ]
        for text, expected in pairs:
            result = html.sanitize(text, user=None)
            TRUE(expected in result)

    def test_smart_mode(self):
        """
        Testing SMART_MODE notifications.
        """
        # Set settings to SMART_MODE
        self.jane.profile.message_prefs = settings.SMART_MODE
        self.jane.profile.save()

        start1 = snapshot()
        # Jane creates a post and a reply.
        post = create_post(user=self.jane)
        reply = create_post(user=self.jane, parent=post)
        end1 = snapshot()

        # Jane has an EMAIL subscription but in SMART mode no emails are
        # for posts that the users creates.
        self.TRUE(PostSub.objects.filter(user=self.jane, post=post, type=settings.EMAIL_TRACKER))
        self.EQ(start1.sub_count + 1, end1.sub_count)
        self.EQ(start1.email_count, end1.email_count)
        self.EQ(start1.msg_count, end1.msg_count)

    def test_email_tracker(self):
        """
        Testing EMAIL_TRACKER notifications.
        """
        # Set settings to SMART_MODE
        self.jane.profile.message_prefs = settings.EMAIL_TRACKER
        self.jane.profile.save()

        start1 = snapshot()
        # Jane creates a post and a reply.
        post = create_post(user=self.jane)
        reply = create_post(user=self.jane, parent=post)
        end1 = snapshot()

        # When EMAIL__TRACKER is enabled a user receives email on their own posts as well.
        self.TRUE(PostSub.objects.filter(user=self.jane, post=post, type=settings.EMAIL_TRACKER))
        self.EQ(start1.sub_count + 1, end1.sub_count)
        self.EQ(start1.email_count + 2, end1.email_count)
        self.EQ(start1.msg_count, end1.msg_count)

    def test_email_on_followups(self):
        """
        Testing email and message notifications on getting an answer or comment.
        """
         # Set settings to SMART_MODE
        self.jane.profile.message_prefs = settings.SMART_MODE
        self.jane.profile.save()

        start1 = snapshot()
        # Jane creates a post.
        post = create_post(user=self.jane)
        reply = create_post(user=self.joe, parent=post)
        comment = create_post(user=self.joe, parent=reply)
        end1 = snapshot()

        # When EMAIL__TRACKER is enabled a user receives email on their own posts as well.
        self.TRUE(PostSub.objects.filter(user=self.jane, post=post, type=settings.EMAIL_TRACKER))

        # Both users now have subscriptions.
        self.EQ(start1.sub_count + 2, end1.sub_count)

        # Janes gets emails and messages.
        self.EQ(start1.email_count + 2, end1.email_count)
        self.EQ(start1.msg_count + 2, end1.msg_count)
        self.EQ(Message.objects.filter(user=self.jane).count(), 2)

    def test_user_tagging(self):
        """
        Testing email notification on watched tags.
        """
        # Joe creates two posts where he tags jane
        self.jane.profile.watched_tags = "hello, jane"
        self.jane.profile.set_tags()

        start1 = snapshot()
        post = create_post(user=self.joe)
        answ = create_post(user=self.joe, parent=post, content="@jane ok")
        comm = create_post(user=self.joe, parent=answ, content="Hello @jane!")
        end1 = snapshot()

        # User should get only one subscription.
        self.EQ(PostSub.objects.filter(post=answ.root, user=self.jane).count(), 1)

        # Jane should get emails because the posts are matching her followed tags
        self.EQ(start1.email_count + 2, end1.email_count)