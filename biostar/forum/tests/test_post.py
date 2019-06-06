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
        self.post = auth.create_post(title="Test", author=self.owner, content="Test",
                                     post_type=models.Post.QUESTION)
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
        self.process_response(response=response)

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

        comment = auth.create_post(title="Test", author=self.owner, content="Test",
                                   post_type=models.Post.COMMENT, root=self.post,
                                   parent=self.post)
        comment2 = auth.create_post(title="Test", author=self.owner, content="Test",
                                   post_type=models.Post.COMMENT, root=self.post,
                                   parent=comment)

        url = reverse("post_view", kwargs=dict(uid=self.post.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.post_view(request=request, uid=self.post.uid)

        self.assertTrue(response.status_code == 200, 'Error rendering comments')

    def test_ajax_subs(self):

        for stype in ["unfollow", "messages", "email", "all", "default"]:

            data = {"sub_type": stype, "root_uid": self.post.uid}
            request = fake_request(url=reverse('vote'), data=data, user=self.owner)
            response = ajax.ajax_subs(request)
            self.assertEqual(response.status_code, 200, f"Could not preform subscription action:{stype}.")

    def preform_votes(self, post, user):
        for vtype in ["upvote", "bookmark", "accept"]:

            data = {"vote_type": vtype, "post_uid": post.uid}
            request = fake_request(url=reverse('vote'), data=data, user=user)
            response = ajax.ajax_vote(request)
            self.assertEqual(response.status_code, 200, f"Could not preform vote:{vtype}.")

    def test_ajax_vote(self):
        """Test the ajax voting using POST request """
        # Create a different user to vote with
        user2 = User.objects.create(username="user", email="user@tested.com", password="tested")

        answer = auth.create_post(title="answer", author=user2, content="tested foo bar too for",
                                  post_type=models.Post.ANSWER, parent=self.post)

        self.preform_votes(post=answer, user=self.owner)
        self.preform_votes(post=self.post, user=self.owner)
        self.preform_votes(post=self.post, user=user2)

        return

    def test_edit_post(self):
        """
        Test post edit for root and descendants
        """
        url = reverse("post_edit", kwargs=dict(uid=self.post.uid))

        # Create a child post to test short form edit
        # Create an existing tested post
        child = auth.create_post(title="Test", author=self.owner, content="Test",
                                 post_type=models.Post.COMMENT, parent=self.post)

        title = "Test title for long test"
        tag_val = "foo,bar,foo"
        content = "Test the content with more things "

        longform_data = dict(title=title, tag_val=tag_val, content=content, post_type=models.Post.TUTORIAL)
        shortform_data = dict(content=content, parent_uid=self.post.uid)

        longform_request = fake_request(url=url, data=longform_data, user=self.owner)
        longform_response = views.edit_post(request=longform_request, uid=self.post.uid)
        self.process_response(longform_response)

        shortform_request = fake_request(url=url, data=shortform_data, user=self.owner)
        shortform_response = views.edit_post(request=shortform_request, uid=child.uid)
        self.process_response(shortform_response)

    def test_post_answer(self):
        """
        Test submitting answer through the post view
        """
        url = reverse("post_answer", kwargs=dict(uid=self.post.uid))

        # Get form data
        data = dict(content="testing answer", parent_uid=self.post.uid)
        request = fake_request(url=url, data=data, user=self.owner)
        response = views.new_answer(request=request, uid=self.post.uid)
        self.process_response(response)
        return

    def test_markdown(self):
        "Test the markdown rendering"
        from django.core import management

        management.call_command("test_markdown")

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")



