import logging
import os
from django.test import TestCase
from django.conf import settings
from biostar.forum import models, markdown
from biostar.accounts.models import User

logger = logging.getLogger('engine')

TEST_CASES = [

    # Top level Post anchors
    ("http://localhost:8000/p/1/", '<p><a href="http://localhost:8000/p/1/">Test</a></p>'),

    # Non-toplevel post anchors
    ("http://localhost:8000/p/1/#2", '<p><a href="http://localhost:8000/p/1/#2">Comment: Test</a></p>'),

    # User profile url pattern
    ("http://localhost:8000/accounts/profile/5", '<p><a href="http://localhost:8000/accounts/profile/5">USER: tested2</a></p>'),

    # Twitter link
    ("https://twitter.com/Linux/status/2311234267", '<p><blockquote class="twitter-tweet"><p lang="en" dir="ltr">w00t! 10,000 followers!</p>&mdash; Linux (@Linux) <a href="https://twitter.com/Linux/status/2311234267?ref_src=twsrc%5Etfw">June 24, 2009</a></blockquote><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script></p>'),

    # Youtube link
    #("https://www.youtube.com/watch?v=dQw4w9WgXcQ", '<p><a href="http://localhost:8000/p/1/#2">Comment: Test</a></p>'),

    # Gist link

    # <code> blocks
    ("""```print &><123http://www.psu.edu```""", """<p><code>print &amp;&gt;&lt;123http://www.psu.edu</code></p>"""),

    # Mentioned user
    ("@test", '<p><a href="/recipes/accounts/profile/5/">tested2</a></p>')



    # Test unclosed tags


    # Test tags and attributes that make the

]


class MarkdownTest(TestCase):
    def setUp(self):
        # Create user
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username="test", email="tested2@tested.com", password="tested")

        self.owner.profile.uid = "5"
        self.owner.profile.save()
        self.owner.save()

        # Create posts used for anchors
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                               type=models.Post.QUESTION, uid="1")

        self.answer = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                                 type=models.Post.ANSWER, uid="2")

        pass

    def test_markdown(self):

        for test in TEST_CASES:
            given, expected = test

            html = markdown.parse(given, clean=True, escape=False)
            html = html.replace("\n", "")
            self.assertEqual(html, expected, f"Error with markdown parsing. input={given}, expexted={expected}, html={html}")