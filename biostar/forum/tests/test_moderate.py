import logging

from django.test import TestCase
from django.urls import reverse

from biostar.accounts.models import User
from unittest.mock import patch, MagicMock
from biostar.forum import models, views
from biostar.forum.moderate import *
from biostar.utils.helpers import fake_request
from biostar.forum.util import get_uuid

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com",
                                         password="tested", is_superuser=True, is_staff=True)
        self.user2 = User.objects.create(username=f"test{get_uuid(10)}", email="test@test.com",
                                         password="test", is_superuser=True, is_staff=True)
        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test", type=models.Post.QUESTION, uid='foo')
        self.uid = 'foo'

        self.owner.save()
        self.post.save()

        pass

    def moderate(self, choices, post, extra={}):

        for action in choices:
            data = {"action": action}
            data.update(extra)

            url = reverse('post_moderate', kwargs=dict(uid=post.uid))
            request = fake_request(url=url, data=data, user=self.owner)
            response = post_moderate(request=request, uid=post.uid)
            self.process_response(response)

        return

    def test_toplevel_moderation(self):
        "Test top level post moderation."
        # Test every moderation action
        choices = ['bump', 'open', 'delete', 'close', 'offtopic', 'relocate']

        self.post = models.Post.objects.create(title="Test",
                                               author=self.owner, content="Test",
                                               type=models.Post.QUESTION, uid='bar')

        self.moderate(choices=choices, post=self.post)

    def test_answer_moderation(self):
        "Test answer moderation."
        choices = ['open', 'delete', 'close', 'offtopic', 'relocate']

        # Create an answer to moderate
        answer = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                  type=models.Post.ANSWER, root=self.post, uid='foo2',
                                  parent=self.post)

        # Add the same amount of giving of the
        self.moderate(choices=choices, post=answer)

        return

    def test_comment_moderation(self):
        "Test comment moderation."
        choices = [ 'open', 'delete', 'close', 'offtopic', 'relocate']

        # Create a comment to moderate
        comment = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post, uid='foo3',
                                   parent=self.post)

        self.moderate(choices=choices, post=comment, extra={'pid': self.post.uid})

    def test_merge_profile(self):
        "Test merging two profiles"

        # Create fake request
        data = {'main': self.owner.email, 'alias': self.user2.email}

        request = fake_request(url=reverse('merge_profile'), data=data, user=self.owner)
        response = views.merge_profile(request=request)

        self.process_response(response)

        pass

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after tested :\nresponse:{response}")



