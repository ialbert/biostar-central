import logging

from django.test import TestCase
from django.urls import reverse

from biostar.forum import models, views, auth, forms
from biostar.engine.test.util import fake_request, get_uuid
from biostar.accounts.models import User

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"tested{get_uuid(10)}", email="tested@tested.com", password="tested")

        # Create an existing tested post
        self.post = auth.create_post(title="Test", author=self.owner, content="Test",
                                     post_type=models.Post.QUESTION)

        self.owner.save()
        pass

    def test_moderate(self):

        # Test every moderation action
        choices = dict(forms.PostModForm.CHOICES)
        for action, wording in choices.items():
            data = {"action": action}
            url = reverse('post_moderate', kwargs=dict(uid=self.post.uid))
            request = fake_request(url=url, data=data, user=self.owner)
            response = views.post_moderate(request=request, uid=self.post.uid)
            self.process_response(response)

        return

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect after tested :\nresponse:{response}")



