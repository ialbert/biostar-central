from __future__ import absolute_import, division, print_function, unicode_literals

import logging, random, re

from django.conf import settings
from django.core.urlresolvers import reverse
from django.core import mail

from biostar3.forum.models import *
from biostar3.forum import auth
from biostar3.tests.client_test_case import ClientTestCase

from faker import Factory

ADMIN_NAME, ADMIN_EMAIL = settings.ADMINS[0]

logging.disable(logging.INFO)

faker = Factory.create()


def add_random_content():
    user = random.choice(User.objects.all())
    parent = random.choice(Post.objects.all())
    content = faker.bs()
    post = auth.create_content_post(content=content, parent=parent, user=user, post_type=None)
    return post


class ClientTests(ClientTestCase):
    def create_post(self, user, parent=None, post_type=None):

        post_type = post_type or Post.QUESTION

        title = faker.sentence()[:100]
        content = faker.text()
        tag1, tag2 = faker.word(), faker.word()
        tag_val = "%s, %s" % (tag1, tag2)

        if parent:
            data = dict(content=content)
            if post_type == Post.ANSWER:
                r = self.post("new_answer", kwargs=dict(pk=parent.id), data=data, follow=True, pattern=content)
            else:
                r = self.post("new_comment", kwargs=dict(pk=parent.id), data=data, follow=True, pattern=content)
            post = Post.objects.filter(content=content).first()
        else:
            data = dict(title=title, tag_val=tag_val, content=content, type=post_type, site=settings.SITE_ID)
            r = self.post("new_post", data=data, follow=True, pattern=title)
            post = Post.objects.filter(title=title).first()

        # Post must exists.
        self.assertTrue(post)

        if post.is_toplevel:
            self.assertTrue(post.author == user)
            self.assertTrue(post.root.author == user)

        return post

    def make_user(self):
        email = faker.email()
        data = dict(email=email, password=email, signup="1")
        r = self.post("sign_up", data=data, follow=True, pattern="success")
        user = User.objects.filter(email=email).first()
        self.assertTrue(user)
        return user

    def test_anonymous_navigation(self):
        """
        Tests client navigation.
        """

        urls = "home user_list tag_list unanswered".split()
        for url in urls:
            r = self.get(url)

        r = self.post("search", data={'q': 'blast'}, follow=True)

    def test_user_navigation(self):
        user = self.make_user()

        # Check a few access pages.
        urls = "home me my_bookmarks my_messages".split()
        for url in urls:
            r = self.get(url, follow=True, pattern=user.name)



    def test_content(self):
        """
        Stress test
        """
        EQ = self.assertEqual

        # Sign up some users. Each makes a post.
        for step in range(10):
            user = self.make_user()
            post = self.create_post(user=user)

        for step in range(25):
            post = add_random_content()
            self.get("post_view", kwargs=dict(pk=post.id), follow=True)

    def test_user_posting(self):
        """
        Test user signup
        """
        EQ, TRUE = self.assertEqual, self.assertTrue

        r = self.get("account_login", pattern='login')

        # Sign up a user.
        jane = self.make_user()

        # User gets a welcome email.
        EQ(len(mail.outbox), 1)

        before = after = len(mail.outbox)

        # Create a question.
        question = self.create_post(user=jane, parent=None)

        # No emails should be sent at this point.
        EQ(before, after)

        # No messages should be sent at this point.
        EQ(Message.objects.all().count(), 0)

        # A different user adds content.
        joe = self.make_user()

        # Joe gets an email
        EQ(len(mail.outbox), before + 1)

        before = len(mail.outbox)
        answer = self.create_post(user=joe, parent=question)

        # Jane gets a message.
        EQ(Message.objects.all().count(), 1)
        TRUE(Message.objects.filter(user=jane).first())

    def test_awards(self):
        """
        Tests user award creation.
        """
        from biostar3.forum import awards

        # Override constants to allow testing them.
        awards.CENTURION_COUNT = 1
        awards.ORACLE_COUNT = 1
        awards.VOTER_COUNT = 1
        awards.SUPPORTER_COUNT = 1

        r = self.get("account_login", pattern='login')

        jane = self.make_user()
        jane.score = awards.CYLON_COUNT
        jane.save()

        jane.profile.info = faker.text()
        jane.profile.save()

        # Create an upvoted question
        question = self.create_post(user=jane, parent=None)
        question.vote_count = awards.CYLON_COUNT
        question.view_count = awards.EPIC_COUNT
        question.subs_count = awards.PROPHET_COUNT

        question.save()

        # Create an upvoted answer.
        answ = self.create_post(user=jane, parent=question, post_type=Post.ANSWER)
        answ.vote_count = awards.CYLON_COUNT
        answ.has_accepted = True
        answ.book_count = awards.BOOKMARK_COUNT
        answ.save()

        # Create an upvoted comment.
        comm = self.create_post(user=jane, parent=question, post_type=Post.COMMENT)
        comm.vote_count = awards.CYLON_COUNT
        comm.reply_count = awards.PUNDIT_COUNT
        comm.save()

        john = self.make_user()
        data = dict(post_id=question.id, vote_type="upvote")
        self.post("vote_submit", data=data, follow=True)

        # The vote has been created.
        self.assertTrue(
            Vote.objects.filter(author=john, type=Vote.UP, post=question).first()
        )

        # Set it up to award each once.
        all_awards = awards.get_awards()

        message_before = Message.objects.all().count()

        for award in all_awards:
            award.check()

        # Each award generated once.
        award_count = Award.objects.all().count()
        message_count = Message.objects.all().count()

        if award_count != len(all_awards):
            # See what is missing.
            for award in Award.objects.all():
                print(award.uuid)

        # User gets one message for each award
        self.assertEqual(award_count, len(all_awards))
        self.assertEqual(message_count, message_before + award_count)

        # Award run the second time around.
        for award in all_awards:
            award.check()

        # No awards should be generated.
        award_count = Award.objects.all().count()
        self.assertEqual(award_count, len(all_awards))







