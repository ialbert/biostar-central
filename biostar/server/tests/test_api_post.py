import json
import logging

from django.test import TestCase
from django.core.urlresolvers import reverse

from ..api import datetime_to_iso
from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User


logging.disable(logging.WARNING)
haystack_logger = logging.getLogger('haystack')


class ApiPostTest(TestCase):
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
        self.post = Post(title=title, content=content, tag_val=tag_val, author=self.user,
                         type=post_type)
        self.post.save()

        # Create a vote.
        self.vote = Vote.objects.create(author=self.user, post=self.post, type=Vote.UP)

    def test_invalid_post_id(self):
        r = self.client.get(reverse('api-post', kwargs={'id': 15}))
        self.assertEqual(r.status_code, 404)

    def test_naked_post(self):
        """
        A naked post with no replies, no bookmarks, no comments, ...
        """
        r = self.client.get(reverse('api-post', kwargs={'id': self.post.id}))
        content = json.loads(r.content)

        expected_post = {
            "answer_count": self.post.root.reply_count,
            "author": "test",
            "author_id": self.user.id,
            "book_count": self.post.book_count,
            "comment_count": self.post.comment_count,
            "creation_date": datetime_to_iso(self.post.creation_date),
            "has_accepted": self.post.has_accepted,
            "id": self.post.id,
            "lastedit_date": datetime_to_iso(self.post.lastedit_date),
            "lastedit_user_id": self.user.id,
            "parent_id": self.post.id,
            "rank": float(self.post.rank),
            "reply_count": self.post.reply_count,
            "root_id": self.post.id,
            "status": self.post.get_status_display(),
            "status_id": self.post.status,
            "subs_count": 1,
            "tag_val": "tag_val",
            "thread_score": self.post.thread_score,
            "title": self.post.title,
            "type": self.post.get_type_display(),
            "type_id": self.post.type,
            "url": 'http://example.com{}'.format(self.post.get_absolute_url()),
            "view_count": self.post.view_count,
            "vote_count": self.post.vote_count,
            "xhtml": self.post.content
        }
        self.assertDictEqual(content, expected_post)

    def test_post(self):
        """
        Regular post with replies, bookmarks, comments, ...
        """
        self.post.reply_count = 3
        self.post.book_count = 4
        self.post.comment_count = 5
        self.post.subs_count = 6
        self.post.view_count = 9
        self.post.vote_count = 7
        self.post.thread_score = 8
        self.post.has_accepted = True
        self.post.rank = 5.5
        self.post.save()

        expected_post = {
            "answer_count": self.post.root.reply_count,
            "author": "test",
            "author_id": self.user.id,
            "book_count": self.post.book_count,
            "comment_count": self.post.comment_count,
            "creation_date": datetime_to_iso(self.post.creation_date),
            "has_accepted": self.post.has_accepted,
            "id": self.post.id,
            "lastedit_date": datetime_to_iso(self.post.lastedit_date),
            "lastedit_user_id": self.user.id,
            "parent_id": self.post.id,
            "rank": float(self.post.rank),
            "reply_count": self.post.reply_count,
            "root_id": self.post.id,
            "status": self.post.get_status_display(),
            "status_id": self.post.status,
            "subs_count": self.post.subs_count,
            "tag_val": "tag_val",
            "thread_score": self.post.thread_score,
            "title": self.post.title,
            "type": self.post.get_type_display(),
            "type_id": self.post.type,
            "url": 'http://example.com{}'.format(self.post.get_absolute_url()),
            "view_count": self.post.view_count,
            "vote_count": self.post.vote_count,
            "xhtml": self.post.content
        }

        r = self.client.get(reverse('api-post', kwargs={'id': self.post.id}))
        content = json.loads(r.content)

        self.assertDictEqual(content, expected_post)