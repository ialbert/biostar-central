import logging
import json
from django.test import TestCase
from django.urls import reverse
from unittest.mock import patch, MagicMock
from biostar.accounts.models import User, Profile

from biostar.forum import models, views, auth, forms, const, ajax
from biostar.utils.helpers import fake_request
from biostar.forum.util import get_uuid

logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested",
                                         is_superuser=True)

        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                               type=models.Post.QUESTION)
        self.owner.save()
        pass

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

        answer = models.Post.objects.create(title="answer", author=user2, content="tested foo bar too for",
                                  type=models.Post.ANSWER, parent=self.post)

        self.preform_votes(post=answer, user=self.owner)
        self.preform_votes(post=self.post, user=self.owner)
        self.preform_votes(post=self.post, user=user2)

    def test_drag_and_drop(self):
        """
        Test AJAX function used to drag and drop.
        """

        # Create comment to self.parent
        comment1 = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post,
                                   parent=self.post)
        # Create another post with parent being above comment
        comment2 = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                   type=models.Post.COMMENT, root=self.post,
                                   parent=comment1)
        comment3 = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                              type=models.Post.COMMENT, root=self.post,
                                              parent=comment2)

        # Move comment3 to comment1
        data = {'parent': comment1.uid, 'uid': comment3.uid}
        url = reverse('drag_and_drop')
        request = fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.drag_and_drop(request)
        self.process_response(json_response)

    def test_digest(self):
        """
        Test AJAX function that toggles users digest options
        """

        data = {'pref': 'daily'}

        url = reverse('ajax_digest')
        request = fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.ajax_digest(request)
        self.process_response(json_response)

    @patch('biostar.forum.models.Post.save', MagicMock(name="save"))
    def test_inplace_edit(self):
        """
        Test AJAX function for inplace edits
        """
        data = {'content': 'New content here foo bar bar', 'title': 'New title here foo bar bar',
                'post_type': models.Post.FORUM, 'tag_val': ['new', 'tag']}

        url = reverse('ajax_edit', kwargs=dict(uid=self.post.uid))

        request = fake_request(url=url, data=data, user=self.owner)
        json_response = ajax.ajax_edit(request, uid=self.post.uid)
        response_data = json.loads(json_response.content)

        self.process_response(json_response)

    def test_disable(self):
        """
        Test project disabling function
        """

        url = reverse('email_disable', kwargs=dict(uid=self.owner.pk))

        request = fake_request(url=url, data={}, user=self.owner)
        json_response = ajax.email_disable(request, uid=self.owner.pk)
        self.process_response(json_response)

    def test_inplace_create(self):
        """
        Test AJAX function used to preform
        """

        nontoplevel_data = {'content': "Addin a new comment foo bar",
                            'parent': self.post.uid}

        url = reverse('ajax_comment_create')

        request = fake_request(url=url, data=nontoplevel_data, user=self.owner)
        nontoplevel_response = ajax.ajax_comment_create(request)
        self.process_response(nontoplevel_response)

    def test_inplace_form(self):
        """
        Test AJAX function used to render inplace form.
        """
        data = {'parent': self.post.uid, 'uid': self.post.uid, 'top': '1'}
        url = reverse('inplace_form')

        request = fake_request(url=url, data=data, user=self.owner, method='GET')
        toplevel_response = ajax.inplace_form(request)
        self.process_response(toplevel_response)

    def test_similar_posts(self):
        """
        Test AJAX function used to produce the similar posts sidebar.
        """
        data = {'query': 'Test'}
        url = reverse('similar_posts', kwargs=dict(uid=self.post.uid))

        request = fake_request(url=url, data=data, user=self.owner, method='GET')
        toplevel_response = ajax.similar_posts(request, uid=self.post.uid)
        self.process_response(toplevel_response)

    def process_response(self, response):
        "Check the response on POST request is redirected"

        response_data = response.content
        response_data = json.loads(response_data)

        self.assertEqual(response_data['status'], 'success', f'Error:{response_data["msg"]}')



