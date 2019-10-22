import logging
from django.test import TestCase
from django.urls import reverse
from biostar.forum import models, views, auth, ajax
from biostar.forum.tests.util import fake_request
from biostar.accounts.models import User

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested")

        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                     type=models.Post.QUESTION)
        self.owner.save()
        pass

    def test_post_create(self):
        """Test post creation with POST request"""

        content = f"@{self.owner.username} " + "testing " * 10

        # Create fake request
        data = {'post_type': models.Post.QUESTION,
                'title': 'title tested post',
                "tag_val": "tested,test3",
                "content": content
                }

        request = fake_request(url=reverse('post_create'), data=data, user=self.owner)
        response = views.new_post(request=request)
        #self.process_response(response=response)

    def test_comment(self):
        """Test adding comment using POST request"""

        data = {"parent_uid": self.post.uid, "content": "tested content for a question"}
        url = reverse('create_comment', kwargs=dict(uid=self.post.uid))

        request = fake_request(url=url, data=data, user=self.owner)
        response = views.new_comment(request=request, uid=self.post.uid)

        self.assertEqual(response.status_code, 302, f"Could not add comments")

    def test_comment_traversal(self):
        """Test comment rendering pages"""

        # Create a couple of comments to traverse

        comment = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post,
                                   parent=self.post)
        comment2 = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post,
                                   parent=comment)

        url = reverse("post_view", kwargs=dict(uid=self.post.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.post_view(request=request, uid=self.post.uid)

        self.assertTrue(response.status_code == 200, 'Error rendering comments')

    def test_edit_post(self):
        """
        Test post edit for root and descendants
        """
        url = reverse("post_edit", kwargs=dict(uid=self.post.uid))

        title = "Test title for long test"
        tag_val = "foo,bar,foo"
        content = "Test the content with more things "

        longform_data = dict(title=title, tag_val=tag_val, content=content, post_type=models.Post.TUTORIAL)

        longform_request = fake_request(url=url, data=longform_data, user=self.owner)
        longform_response = views.edit_post(request=longform_request, uid=self.post.uid)
        #self.process_response(longform_response)

    def test_post_answer(self):
        """
        Test submitting answer through the post view
        """
        url = reverse("post_view", kwargs=dict(uid=self.post.uid))

        # Get form data
        data = dict(content="testing answer", parent_uid=self.post.uid)
        request = fake_request(url=url, data=data, user=self.owner)
        response = views.post_view(request=request, uid=self.post.uid)
        #self.process_response(response)
        return

    def test_markdown(self):
        "Test the markdown rendering"
        from django.core import management

        management.call_command("test_markdown")

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")



