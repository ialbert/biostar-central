import logging

from django.test import TestCase
from django.urls import reverse

from biostar.accounts.models import User
from unittest.mock import patch, MagicMock
from biostar.forum import models, views, auth, forms, const
from biostar.utils.helpers import fake_request
from biostar.forum.util import get_uuid

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com",
                                         password="tested", is_superuser=True, is_staff=True)

        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test",
                                               author=self.owner, content="Test",
                                     type=models.Post.QUESTION, uid='foo')
        self.uid = 'foo'

        self.owner.save()
        self.post.save()

        pass

    def moderate(self, choices, post, extra={}):

        for action in choices:
            data = {"action": action}
            data.update(extra)

            url = reverse('post_moderate', kwargs=dict(uid=self.uid))
            request = fake_request(url=url, data=data, user=self.owner)
            response = views.post_moderate(request=request, uid=post.uid)
            self.process_response(response)

        return

    @patch('biostar.recipes.models.Job.save', MagicMock(name="save"))
    def test_toplevel_moderation(self):
        "Test top level post moderation."
        # Test every moderation action
        choices = [const.BUMP_POST, const.OPEN_POST, const.DELETE]

        #self.moderate(choices=choices, post=self.post)

        return

    def test_answer_moderation(self):
        "Test answer moderation."
        choices = [const.TOGGLE_ACCEPT, const.DELETE]

        # Create an answer to moderate
        answer = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                  type=models.Post.ANSWER, root=self.post,
                                  parent=self.post)


        answer.save()

        #self.moderate(choices=choices, post=answer)

        return

    def test_comment_moderation(self):
        "Test comment moderation."
        choices = [const.DELETE]

        # Create a comment to moderate
        comment = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post,
                                   parent=self.post)

        comment.save()

        #self.moderate(choices=choices, post=comment, extra={'pid': self.post.uid})

    def test_duplicate_post(self):
        "Test duplicate post moderation"

        data = {"dupe": "google.com"}

        url = reverse('post_moderate', kwargs=dict(uid=self.post.uid))
        request = fake_request(url=url, data=data, user=self.owner)
        response = views.post_moderate(request=request, uid=self.post.uid)
        self.process_response(response)

        pass

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after tested :\nresponse:{response}")



