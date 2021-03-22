import logging
import os
import shutil
from django.core import management
from django.urls import reverse
from django.test import TestCase, override_settings
from django.conf import settings
from biostar.forum import models, views, search, tasks, feed
from biostar.utils.helpers import fake_request
from biostar.accounts.models import User

logger = logging.getLogger('engine')

TEST_DATABASE_NAME = f"test_{settings.DATABASE_NAME}"
TEST_DEBUG = True

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'test'))
TEST_INDEX_DIR = TEST_ROOT
TEST_INDEX_NAME = "index"


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested")
        self.staff_user = User.objects.create(username=f"test2", is_superuser=True, is_staff=True,
                                              email="tested@staff.com", password="tested")

        # Create an existing tested post
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                     type=models.Post.QUESTION)
        self.owner.save()
        pass

    @override_settings(SEND_MAIL=True)
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

    def test_user_create_task(self):
        """
        Test task used to create user awards
        """

        self.owner.profile.text = "TESTING" * 1000
        self.owner.profile.score = 1000
        self.owner.profile.save()
        tasks.create_user_awards(self.owner.id)


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

    def Xtest_edit_post(self):
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

    @override_settings(SEND_MAIL=True)
    def test_post_answer(self):
        """
        Test submitting answer through the post view
        """
        url = reverse("post_view", kwargs=dict(uid=self.post.uid))

        # Get form data
        data = dict(content="testing answer", parent_uid=self.post.uid)
        request = fake_request(url=url, data=data, user=self.owner)
        response = views.post_view(request=request, uid=self.post.uid)
        return

    @override_settings(DEBUG=TEST_DEBUG)
    def test_populate(self):
        "Test forum populating "

        management.call_command('populate', n_users=10, n_messages=10, n_votes=10, n_posts=10)

    def test_markdown(self):
        "Test the markdown rendering"
        from django.core import management

        #management.call_command("test_markdown")

    def process_response(self, response):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect :\nresponse:{response}")


@override_settings(INDEX_DIR=TEST_INDEX_DIR, INDEX_NAME=TEST_INDEX_NAME, DATABASE_NAME=TEST_DATABASE_NAME)
class PostSearchTest(TestCase):

    def setUp(self):
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested")

        # Delete test search index on each start up.
        if os.path.exists(TEST_INDEX_DIR):
            shutil.rmtree(TEST_INDEX_DIR)

        # Create some posts to index.
        self.limit = 10
        for p in range(self.limit):
            # Create an existing tested post
            self.post = models.Post.objects.create(title=f"Test post-{p} ", author=self.owner,
                                                   content=f"Test post-{p} ",
                                                   type=models.Post.QUESTION)
        self.owner.save()

        # Crawl through posts and create test index.
        search.crawl(reindex=True, overwrite=True, limit=1000)

        # There should not be any more unidexed posts.
        unindexed_posts = models.Post.objects.filter(indexed=False).exists()
        #TODO: put back in
        #self.assertFalse(unindexed_posts, "Posts not correctly indexed.")

    def test_search_view(self):
        """
        Test search view
        """

        url = reverse("post_search")

        # Get form data
        data = dict(query="Test post")
        request = fake_request(url=url, data=data, method="GET", user=self.owner)

        # Preform search on test index
        response = views.post_search(request=request)

        self.assertEqual(response.status_code, 200, "Error preforming search.")

    def test_search_funcs(self):
        """
        Test functions associated with search.
        """
        query = "Test"
        whoosh_search = search.perform_search(query)

        search.print_info()
        # TODO: put back in
        #self.assertTrue(len(whoosh_search), f"Whoosh search returned no results. At least {self.limit} expected")