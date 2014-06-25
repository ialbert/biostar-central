import logging

from django.test import TestCase
from django.core.urlresolvers import reverse

from biostar.apps.posts.models import Post, Tag
from biostar.apps.users.models import User


logging.disable(logging.WARNING)
haystack_logger = logging.getLogger('haystack')


class ApiPostTest(TestCase):
    def setUp(self):
        # Disable haystack logger (testing will raise errors on more_like_this field in templates).
        haystack_logger.setLevel(logging.CRITICAL)

        # Create a user.
        self.user = User.objects.create(email='test@test.com', password='...')

        # Create tags.
        t1 = Tag.objects.create(name='mytag1')
        t2 = Tag.objects.create(name='mytag2')

        # Create a post.
        title = "Post 1, title needs to be sufficiently long"
        content = ('Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod'
                   'tempor incididunt ut labore et dolore magna aliqua.')
        post_type = Post.QUESTION
        self.post = Post(title=title, content=content, author=self.user, type=post_type)
        self.post.save()
        self.post.tag_set.add(t1)
        self.post.tag_set.add(t2)
        self.post.save()

    def test_home_page(self):
        """
        Ensure that the new post in the home page uses the right tags and links.
        """
        r = self.client.get(reverse('home'))
        self.assertIn('<a class="tag" href="/t/mytag1/">mytag1</a>', r.content)
        self.assertIn('<a class="tag" href="/t/mytag2/">mytag2</a>', r.content)

    def test_topic_page_single_tag(self):
        """
        Ensure that the new post in the topic page uses the right tags and links when only 1 tag
        is selected.
        """
        r = self.client.get(reverse('topic-list', kwargs={'topic': 'mytag1'}))
        self.assertIn('<a class="tag" href="/t/mytag1/">mytag1</a>', r.content)
        self.assertIn('<a class="tag" href="/t/mytag1+mytag2/">mytag2</a>', r.content)

    def test_topic_page_multiple_tags(self):
        """
        Ensure that the new post in the topic page uses the right tags and links when multiple
        tags are selected.
        """
        r = self.client.get(reverse('topic-list', kwargs={'topic': 'mytag1+mytag2'}))
        self.assertIn('<a class="tag" href="/t/mytag1+mytag2/">mytag1</a>', r.content)
        self.assertIn('<a class="tag" href="/t/mytag1+mytag2/">mytag2</a>', r.content)
        self.assertIn('Filtering by tags: mytag1 OR mytag2', r.content)