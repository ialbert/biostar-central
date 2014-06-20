from datetime import datetime, timedelta
import json
import logging

from django.test import TestCase
from django.core.urlresolvers import reverse

from ..api import datetime_to_iso
from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User


logging.disable(logging.WARNING)
haystack_logger = logging.getLogger('haystack')


class ApiUserTest(TestCase):
    def setUp(self):
        # Disable haystack logger (testing will raise errors on more_like_this field in templates).
        haystack_logger.setLevel(logging.CRITICAL)

        # Create a user.
        self.user = User.objects.create(email='test@test.com', password='...')

        # Edit date_joined.
        self.user.profile.date_joined = datetime.now() - timedelta(days=12)
        self.user.profile.save()

    def test_invalid_user_id(self):
        r = self.client.get(reverse('api-user', kwargs={'id': 15}))
        self.assertEqual(r.status_code, 404)

    def test_user(self):
        r = self.client.get(reverse('api-user', kwargs={'id': self.user.id}))
        content = json.loads(r.content)

        self.assertEqual(content['date_joined'][:10].encode(),
                         datetime_to_iso(self.user.profile.date_joined)[:10])
        self.assertEqual(content['last_login'][:10].encode(),
                         datetime_to_iso(datetime.today())[:10])
        self.assertEqual(content['id'], self.user.id)
        self.assertEqual(content['joined_days_ago'], 12)
        self.assertEqual(content['name'], 'test')
        self.assertEqual(content['vote_count'], 0)

    def test_user_w_vote(self):
        self.create_vote()

        r = self.client.get(reverse('api-user', kwargs={'id': self.user.id}))
        content = json.loads(r.content)
        self.assertEqual(content['vote_count'], 1)

    def create_vote(self):
        post = self.create_a_post()
        Vote.objects.create(author=self.user, post=post, type=Vote.UP)

    def create_a_post(self):
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        post_type = Post.QUESTION
        tag_val = 'tag_val'

        post = Post(title=title, content=content, tag_val=tag_val, author=self.user,
                    type=post_type)
        post.save()

        # Triggers a new post save.
        post.add_tags(post.tag_val)

        return post