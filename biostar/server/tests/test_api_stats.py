import logging
from datetime import datetime, timedelta
import json

from django.test import TestCase
from django.core.urlresolvers import reverse
from django.test.utils import override_settings

from biostar.apps.posts.models import Post, Vote
from biostar.apps.users.models import User


haystack_logger = logging.getLogger('haystack')
# Disable haystack logger (testing will raise errors on more_like_this field in templates).
haystack_logger.setLevel(logging.CRITICAL)
logging.disable(logging.CRITICAL)


# Note: all test here use settings.DEBUG = True to avoid stats files to be created.


class ApiStatsTest1(TestCase):
    @override_settings(DEBUG=True)
    def test_day_no_posts(self):
        """
        No posts in the db.
        """
        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 1}))
        self.assertEqual(r.content, '{}')

    @override_settings(DEBUG=True)
    def test_date_no_posts(self):
        """
        No posts in the db.
        """
        r = self.client.get(reverse('api-stats-on-date',
                                    kwargs={'year': '1900', 'month': '01', 'day': '01'}))
        content = json.loads(r.content)
        expected_data = {
            "answers": 0,
            "comments": 0,
            "date": "1900-01-01T00:00:00",
            "new_posts": [],
            "new_users": [],
            "new_votes": [],
            "questions": 0,
            "timestamp": -2208988800,
            "toplevel": 0,
            "users": 0,
            "votes": 0
        }
        self.assertDictEqual(content, expected_data)


class ApiStatsTest2(TestCase):
    def setUp(self):
        # Create a user and edit the date joined.
        self.user = User.objects.create(email='test@test.com', password='...')
        self.user.profile.date_joined = datetime.today() - timedelta(days=3)
        self.user.profile.save()

        # Create a question and a vote.
        self.post = self.create_post(self.user, Post.QUESTION, days=3)
        Vote.objects.create(author=self.user, post=self.post, type=Vote.UP)

    def create_post(self, user, post_type, days=3):
        # Create a post.
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        tag_val = 'tag_val'
        post = Post(title=title, content=content, tag_val=tag_val, author=user, type=post_type, )
        post.save()
        post.creation_date = datetime.today() - timedelta(days=days)
        post.save()
        return post

    @override_settings(DEBUG=True)
    def test_day_a_post_a_user(self):
        """
        A new post and a new user.
        """
        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 0}))
        content = json.loads(r.content)
        expected_data = {
            "answers": 0,
            "comments": 0,
            "date": content['date'],  # Hard to test cause timezones are involved.
            "new_posts": [
                self.post.id
            ],
            "new_users": [
                self.user.id
            ],
            "new_votes": [],
            "questions": 1,
            "timestamp": content['timestamp'],  # Hard to test cause timezones are involved.
            "toplevel": 1,
            "users": 1,
            "votes": 0
        }
        self.assertDictEqual(content, expected_data)
        # Use the following lines to debug the content of the dictionaries.
        #for key, val in expected_data.items():
        #    print(key, content[key], val)
        #    self.assertEqual(content[key], val)

    @override_settings(DEBUG=True)
    def test_date_a_post_a_user(self):
        """
        A new post and a new user.
        """
        date = datetime.today() - timedelta(days=2)
        r = self.client.get(reverse('api-stats-on-date',
                                    kwargs={
                                        'year': date.year,
                                        'month': str(date.month).zfill(2),
                                        'day': str(date.day).zfill(2)
                                    }))
        content_date = json.loads(r.content)

        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 1}))
        content_day = json.loads(r.content)

        self.assertEqual(content_day, content_date)


class ApiStatsTest3(TestCase):
    def setUp(self):
        # Create a user.
        self.user = User.objects.create(email='test@test.com', password='...')
        self.user.profile.date_joined = datetime.today() - timedelta(days=4)
        self.user.profile.save()

        self.first_question = self.create_post(Post.QUESTION, days=4)
        self.question = self.create_post(Post.QUESTION)
        self.answer = self.create_post(Post.ANSWER, self.question)
        self.comment = self.create_post(Post.COMMENT, self.question)

    def create_post(self, post_type, parent=None, days=3):
        # Create a post.
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        tag_val = 'tag_val'
        post = Post(title=title, content=content, tag_val=tag_val, author=self.user,
                    type=post_type, parent=parent)
        #post.save()
        post.creation_date = datetime.today() - timedelta(days=days)
        post.save()
        return post

    @override_settings(DEBUG=True)
    def test_day_question_comment_answer(self):
        """
        There are an answer, a comment and 2 questions.
        """
        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 3}))
        content = json.loads(r.content)

        self.assertEqual(content['answers'], 1)
        self.assertEqual(content['comments'], 1)
        self.assertEqual(content['questions'], 2)

    @override_settings(DEBUG=True)
    def test_date_question_comment_answer(self):
        """
        There are an answer, a comment and 2 questions.
        """
        date = datetime.today() - timedelta(days=1)
        r = self.client.get(reverse('api-stats-on-date',
                                    kwargs={
                                        'year': date.year,
                                        'month': str(date.month).zfill(2),
                                        'day': str(date.day).zfill(2)
                                    }))
        content_date = json.loads(r.content)

        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 3}))
        content_day = json.loads(r.content)

        self.assertEqual(content_day, content_date)