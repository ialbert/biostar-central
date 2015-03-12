from __future__ import absolute_import, division, unicode_literals

import logging, random, re

from django.test import TestCase
from django.conf import settings
from django.core.urlresolvers import reverse
from django.core import mail
from django.test import Client

from biostar3.forum.models import *
from biostar3.forum import auth

from faker import Factory

HOST = "www.lvh.me:8080"

ADMIN_NAME, ADMIN_EMAIL = settings.ADMINS[0]

logging.disable(logging.ERROR)

faker = Factory.create()


def add_random_content():
    user = random.choice(User.objects.all())
    parent = random.choice(Post.objects.all())
    content = faker.bs()
    auth.create_content_post(content=content, parent=parent, user=user, post_type=None)
    post = Post.objects.filter(content=content)
    return user, parent, post


class ClientTests(TestCase):
    def get(self, client, url, code=200, follow=False, pattern=None):
        r = client.get(reverse(url), follow=follow)
        self.assertEqual(r.status_code, code)
        if pattern:
            self.assertTrue(re.search(pattern, r.content, re.IGNORECASE))
        return r

    def post(self, client, url, kwargs={}, data={}, pattern=None, follow=False):
        r = client.post(reverse(url, kwargs=kwargs), data, follow=follow)

        if follow:
            self.assertEqual(r.status_code, 200)
        else:
            self.assertEqual(r.status_code, 302)

        if pattern:
            result = re.search(pattern, r.content, re.IGNORECASE)
            if not result:
                print "Unable to find %s pattern in content." % (pattern)
                self.assertTrue(result)

        return r

    def test_navigation(self):
        """
        Tests client navigation.
        """

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "home")
        r = self.post(c, "search", data={'q': 'blast'}, follow=True)


    def make_post(self, client, user, parent=None, post_type=None):

        post_type = Post.QUESTION

        title = faker.sentence()[:100]
        content = faker.text()
        tag1, tag2 = faker.word(), faker.word()
        tags = "%s, %s" % (tag1, tag2)

        if parent:
            data = dict(content=content)
            if post_type == Post.ANSWER:
                r = self.post(client, "new_answer", kwargs=dict(pk=parent.id), data=data, follow=True, pattern=title)
            else:
                r = self.post(client, "new_comment", kwargs=dict(pk=parent.id), data=data, follow=True, pattern=title)
            post = Post.objects.filter(content=content).first()
        else:
            data = dict(title=title, tags=tags, content=content, type=post_type)
            r = self.post(client, "new_post", data=data, follow=True, pattern=title)
            post = Post.objects.filter(title=title).first()

        # Post must exists.
        self.assertTrue(post)

        if post.is_toplevel:
            self.assertTrue(post.author == user)
            self.assertTrue(post.root.author == user)

        return post

    def make_user(self, client):
        email = faker.email()
        data = dict(email=email, password=email, signup="1")
        r = self.post(client, "sign_up", data=data, follow=True, pattern="success")
        user = User.objects.filter(email=email).first()
        self.assertTrue(user)
        return user


    def test_content(self):
        """
        Stress test
        """
        EQ = self.assertEqual

        c = Client(HTTP_HOST=HOST)

        # Sign up some users. Each makes a post.
        for step in range(10):
            user = self.make_user(c)
            post = self.make_post(c, user=user)

        add_random_content()

    def test_user_posting(self):
        """
        Test user signup
        """
        EQ = self.assertEqual

        c = Client(HTTP_HOST=HOST)
        r = self.get(c, "account_login", pattern='social authentication')

        # Sign up a user.
        user = self.make_user(c)

        # Check a few access pages.
        for patt in "new_post home me".split():
            r = self.get(c, patt, follow=True)

        before = after = len(mail.outbox)

        # Create a question.
        post = self.make_post(client=c, user=user, parent=None)

        # No emails should be sent at this point.
        EQ(before, after)


        EQ(Message.objects.all().count(), 0)

        # Change notification.






