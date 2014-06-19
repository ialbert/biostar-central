from datetime import datetime
import json
import logging

from django.test import TestCase
from django.core.urlresolvers import reverse

from ..api import datetime_to_iso, datetime_to_unix
from biostar.apps.posts.models import Post
from biostar.apps.users.models import User


logging.disable(logging.WARNING)
haystack_logger = logging.getLogger('haystack')


class ApiTrafficTest(TestCase):
    def setUp(self):
        # Disable haystack logger (testing will raise errors on more_like_this field in templates).
        haystack_logger.setLevel(logging.CRITICAL)

    def test_no_traffic(self):
        """
        There is no posts in the db.
        """
        r = self.client.get(reverse('api-traffic'))
        now = datetime.now()
        content = json.loads(r.content)
        self.assertTrue(content['date'].startswith(datetime_to_iso(now)[:11]))
        self.assertLess(datetime_to_unix(now) - content['timestamp'], 100)
        self.assertEqual(content['post_views_last_60_min'], 0)

    def test_one_visit(self):
        """
        1 Post view.
        """
        # Create a user.
        user = self.create_a_user()

        # Create a post.
        post = self.create_a_post(user)

        # Create a post-view.
        self.client.get(reverse('post-details', kwargs={'pk': post.pk}))

        # Check traffic.
        r = self.client.get(reverse('api-traffic'))
        now = datetime.now()
        content = json.loads(r.content)
        self.assertTrue(content['date'].startswith(datetime_to_iso(now)[:11]))
        self.assertLess(datetime_to_unix(now) - content['timestamp'], 100)
        self.assertEqual(content['post_views_last_60_min'], 1)

    def create_a_user(self):
        with self.settings(CAPTCHA=False, TRUST_VOTE_COUNT=0):
            email_address = 'test@test.com'
            self.client.post(reverse("account_signup"),
                             {
                                 'email': email_address,
                                 'password1': 'password',
                                 'password2': 'password',
                                 'follow': True,
                             },)
        return User.objects.get(email=email_address)

    @staticmethod
    def create_a_post(user):
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        post_type = Post.QUESTION
        tag_val = 'tag_val'

        post = Post(title=title, content=content, tag_val=tag_val, author=user, type=post_type)
        post.save()

        # Triggers a new post save.
        post.add_tags(post.tag_val)

        return post