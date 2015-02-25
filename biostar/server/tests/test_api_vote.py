import json
import logging

from django.test import TestCase
from django.core.urlresolvers import reverse

from ..api import datetime_to_iso
from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User


logging.disable(logging.WARNING)
haystack_logger = logging.getLogger('haystack')


class ApiVoteTest(TestCase):
    def setUp(self):
        # Disable haystack logger (testing will raise errors on more_like_this field in templates).
        haystack_logger.setLevel(logging.CRITICAL)

        # Create a user.
        self.user = User.objects.create(email='test@test.com', password='...')

        # Create a post.
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        post_type = Post.QUESTION
        tag_val = 'tag_val'
        post = Post(title=title, content=content, tag_val=tag_val, author=self.user,
                    type=post_type)
        post.save()

        # Create a vote.
        self.vote = Vote.objects.create(author=self.user, post=post, type=Vote.UP)

    def test_invalid_vote_id(self):
        r = self.client.get(reverse('api-vote', kwargs={'id': 15}))
        self.assertEqual(r.status_code, 404)

    def test_vote(self):
        r = self.client.get(reverse('api-vote', kwargs={'id': self.vote.id}))
        content = json.loads(r.content)

        self.assertEqual(content['author'], 'test')
        self.assertEqual(content['author_id'], self.user.id)
        self.assertEqual(content['date'], datetime_to_iso(self.vote.date))
        self.assertEqual(content['id'], self.vote.id)
        self.assertEqual(content['post_id'], self.vote.post.id)
        self.assertEqual(content['type'], self.vote.get_type_display())
        self.assertEqual(content['type_id'], self.vote.type)
