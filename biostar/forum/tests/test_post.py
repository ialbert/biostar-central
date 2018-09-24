import logging

from django.test import TestCase
from django.urls import reverse

from biostar.forum import models, views, auth, forms
from biostar.engine.test.util import fake_request
from biostar.accounts.models import User


logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username="test", email="test@test.com", password="testing")

        # Create an existing test post
        self.post = auth.create_post(title="Test", author=self.owner, content="Test",
                                     post_type=models.Post.QUESTION)

        self.owner.save()
        pass

    def test_post_create(self):
        """Test post creation with POST request"""

        # Create fake request
        data = {'post_type': models.Post.QUESTION,
                'title': 'title testing post',
                "tag_val": "test",
                "content": "test content for a question"
                }

        request = fake_request(url=reverse('post_create'), data=data, user=self.owner)
        response = views.post_create(request=request)
        self.process_response(response=response)

    def test_comment(self):
        """Test adding comment using POST request"""

        data = {"parent_uid": self.post.uid,"content": "test content for a question"}
        request = fake_request(url=reverse('post_comment'), data=data, user=self.owner)
        response = views.comment(request=request)

        wrong_data = {"content": "test content for a question"}
        wrong_request = fake_request(url=reverse('post_comment'), data=wrong_data, user=self.owner)
        wrong_response = views.comment(request=wrong_request)

        self.process_response(response=response)

        self.assertEqual(wrong_response.url, reverse("post_list"))
        self.process_response(response=wrong_response)

    def make_votes(self, post, user):

        for vtype in ["upvote", "bookmark", "accept"]:

            data = {"vote_type": vtype, "post_uid": post.uid}
            request = fake_request(url=reverse('vote'), data=data, user=user)
            response = views.ajax_vote(request=request)

    def test_vote(self):
        """Test the ajax voting using POST request """
        user2 = User.objects.create(username="user", email="user@test.com", password="test")

        answer = auth.create_post(title="answer", author=user2, content="test foo bar too for",
                                  post_type=models.Post.ANSWER, parent=self.post)

        self.make_votes(post=answer, user=self.owner)
        self.make_votes(post=self.post, user=self.owner)
        self.make_votes(post=self.post, user=user2)

        return

    def test_answer(self):

        return

    def test_moderate(self):

        # Test every moderation action
        for action in forms.PostModForm.CHOICES:
            pass

        return

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after testing :\nresponse:{response}")






