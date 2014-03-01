from django.test import TestCase, SimpleTestCase
from django.test import Client
from django.core.urlresolvers import reverse
from django.conf import settings
from biostar.apps.users.models import User
from biostar.apps.posts.models import Post, Tag, PostView, Subscription
from biostar.apps.messages.models import Message, MessageBody

import logging, random

logging.disable(logging.WARNING)

user_count = lambda: User.objects.all().count()
post_count = lambda: Post.objects.all().count()
subs_count = lambda: Subscription.objects.all().count()
msg_count = lambda: Message.objects.all().count()
get_user = lambda x: User.objects.get(email=x)


class UserTest(TestCase):
    # The name of test users
    NAME_1, EMAIL_1, PASSWD_1 = "John Doe", "user1@example.org", "0123567"
    NAME_2, EMAIL_2, PASSWD_2 = "Jane Doe", "user2@example.org", "3456789"

    USER_DATA = [
        (EMAIL_1, PASSWD_1),
        (EMAIL_2, PASSWD_2),
    ]

    # The name of test posts
    TITLE_1 = "Post 1, title needs to be sufficiently long"
    CAT_1, TAG_VAL_1 = Post.QUESTION, "tagA tagB galaXY"

    TITLE_2 = "Post 2, title needs to be sufficiently long"
    CAT_2, TAG_VAL_2 = Post.JOB, "jobA jobB galaxy"

    CONTENT = """
    Lorem ipsum dolor sit amet, consectetur adipisicing elit,
    sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
    """

    POST_DATA = [
        (TITLE_1, CAT_1, TAG_VAL_1),
        (TITLE_2, CAT_2, TAG_VAL_2),
    ]

    def code(self, response, code=200):
        self.assertEqual(response.status_code, code)

    def setUp(self):
        # Sign up then log out each user.
        for email, passwd in self.USER_DATA:
            self.sign_up(email, passwd)

    def sign_up(self, email, passwd):
        count = user_count()
        r = self.client.post(reverse("account_signup"), dict(email=email, password1=passwd, password2=passwd),
                             follow=True)
        self.assertContains(r, "My Tags")
        self.assertEqual(user_count(), count + 1)
        self.logout()

    def login(self, email, passwd):
        "Logs in a user"
        r = self.client.post(reverse("account_login"), dict(login=email, password=passwd), follow=True)
        self.assertContains(r, "My Tags")


    def logout(self):
        "Logs out the current user."
        r = self.client.get(reverse("logout"), follow=True)
        self.assertNotContains(r, "My Tags")


    def test_user_login(self):
        "Test that each user can log in."
        eq = self.assertEqual
        for email, passwd in self.USER_DATA:
            self.login(email, passwd)
            self.logout()


    def create_new_post(self, title, post_type, tag_val):
        p_count = post_count()
        s_count = subs_count()
        r = self.client.post(
            reverse("new-post"),
            dict(title=title, tag_val=tag_val, post_type=post_type, content=self.CONTENT),
        )

        # Needs to redirect to post
        self.code(r, 302)

        # After creating a new post the post count and subscription counts increase.
        self.assertEqual(post_count(), p_count + 1)
        self.assertEqual(subs_count(), s_count + 1)



    def create_new_answer(self, post):
        r = self.client.post(
            reverse("new-answer", kwargs=dict(pid=post.id)),
            dict(content=self.CONTENT),
        )
        self.code(r, 302)


    def get_post(self, pk):
        "Gets a post and returns it"
        post = Post.objects.get(pk=pk)
        r = self.client.get(reverse("post-details", kwargs=dict(pk=pk)))
        if post.is_toplevel:
            self.assertContains(r, post.title)
        self.code(r)

        # Verify that a subscription exists for this post and author.

        self.assertTrue(Subscription.objects.get_subs(post).filter(user=post.author).count() == 1)
        return post


    def test_user_new_post(self):
        "Test that each user can create a new post."
        eq = self.assertEqual

        for email, passwd in self.USER_DATA:
            self.login(email, passwd)
            for title, post_type, tag_val in self.POST_DATA:
                # Create unique titles
                title = title + email
                self.create_new_post(title=title, post_type=post_type, tag_val=tag_val)
                post = Post.objects.get(title=title)
                self.get_post(post.id)
            self.logout()


    def test_user_answer(self):
        "Test posting an answer."
        self.login(self.EMAIL_1, self.PASSWD_1)
        title = self.TITLE_1
        self.create_new_post(title=title, post_type=self.CAT_1, tag_val=self.TAG_VAL_1)
        post1 = Post.objects.get(title=title)
        post2 = self.get_post(post1.id)
        self.assertEqual(post1, post2)

        # Same user adds a new answer.
        p_count = post_count()
        m_count = msg_count()
        self.create_new_answer(post1)

        # Post counts increased.
        self.assertEqual(post_count(), p_count + 1)

        # No message has been added because it is the same user.
        self.assertEqual(msg_count(), m_count)
        self.logout()

        # A different user adds an answer
        self.login(self.EMAIL_2, self.PASSWD_2)
        self.create_new_answer(post1)

        # Post count increases
        self.assertEqual(post_count(), p_count + 2)

        # A message is added for the author of the parent.
        self.assertEqual(msg_count(), m_count + 1)
        self.assertTrue(Message.objects.filter(user__email=self.EMAIL_1).count() == 1)

        # Test voting and that it applies to user and posts
        user1 = get_user(self.EMAIL_1)
        post = Post.objects.get(title=self.TITLE_1, type=Post.QUESTION)

        # First access adds a vote.
        r = self.client.post(reverse("vote-submit"), data=dict(post_id=post.id, vote_type="vote"))
        user2 = get_user(self.EMAIL_1)
        self.assertEqual(user1.score + 1, user2.score)

        # Seconds access removes a vote.
        r = self.client.post(reverse("vote-submit"), data=dict(post_id=post.id, vote_type="vote"))
        user3 = get_user(self.EMAIL_1)
        self.assertEqual(user1.score, user3.score)

        # Bookmarks also add reputation.
        r = self.client.post(reverse("vote-submit"), data=dict(post_id=post.id, vote_type="bookmark"))
        user4 = get_user(self.EMAIL_1)
        self.assertEqual(user1.score + 1, user4.score)


    def test_post_rendering(self):
        "Stress test."
        emails = [ "%s@test.org" % x for x in range(10) ]
        passwd = "1234567"
        for email in emails:
            self.sign_up(email, passwd)

        top_types = [Post.QUESTION, Post.JOB, Post.FORUM, Post.PAGE]
        for count in range(5):
            email = random.choice(emails)
            post_type = random.choice(top_types)
            self.login(email, passwd)
            self.create_new_post(title=self.TITLE_1, post_type=post_type, tag_val=self.TAG_VAL_1)
            self.logout()

        for count in range(10):
            email = random.choice(emails)
            self.login(email, passwd)
            id = random.choice(range(1, Post.objects.all().count()))
            self.get_post(pk=id)
            post = Post.objects.get(pk=id)
            self.create_new_answer(post)
            self.logout()

class SiteTest(TestCase):

    def code(self, response, code=200):
        self.assertEqual(response.status_code, code)

    def test_site_navigation(self):
        "Testing site navigation."

        eq = self.assertEqual

        # Main site navigation.
        names = "home user-list tag-list help about faq policy rss latest-feed".split()
        for name in names:
            r = self.client.get(reverse(name))
            self.code(r)

        # Check that default categories work.
        for topic in settings.CATEGORIES:
            r = self.client.get(reverse("topic-list", kwargs=dict(topic=topic)))
            self.code(r)

    def test_redirects(self):
        "Testing page redirects."

        # Pages with redirects.
        names = "login logout signup new-post".split()
        for name in names:
            r = self.client.get(reverse(name))
            self.code(r, 302)

        # Pages that take parameters and redirect.
        names = "user-edit post-edit".split()
        for name in names:
            r = self.client.get(reverse(name, kwargs=dict(pk=1)))
            self.code(r, 302)



