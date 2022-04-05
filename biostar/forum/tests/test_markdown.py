import logging
import os
from django.test import TestCase
from django.conf import settings
from biostar.forum import models, markdown
from biostar.accounts.models import User

logger = logging.getLogger('engine')

PORT = ':' + settings.HTTP_PORT if settings.HTTP_PORT else ''
SITE_URL = f"{settings.SITE_DOMAIN}{PORT}"



TEST_CASES = [

    # Top level Post anchors
    (f"{settings.PROTOCOL}://{SITE_URL}/p/1/ ", f'<p><a href="{settings.PROTOCOL}://{SITE_URL}/p/1/" rel="nofollow">Test</a></p>'),

    # Non-toplevel post anchors
    (f"{settings.PROTOCOL}://{SITE_URL}/p/1/#2 ", f'<p><a href="{settings.PROTOCOL}://{SITE_URL}/p/1/#2" rel="nofollow">Comment: Test</a></p>'),

    # User profile url pattern
    (f"{settings.PROTOCOL}://{SITE_URL}/u/5 ", f'<p><a href="{settings.PROTOCOL}://{SITE_URL}/u/5" rel="nofollow">tested2</a></p>'),

    # Twitter link
    ("https://twitter.com/Linux/status/2311234267", '<p><blockquote class="twitter-tweet"><p lang="en" dir="ltr">w00t! 10,000 followers!</p>&mdash; Linux (@Linux) <a href="https://twitter.com/Linux/status/2311234267?ref_src=twsrc%5Etfw">June 24, 2009</a></blockquote><script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script></p>'),

    # Youtube link
    ("https://www.youtube.com/watch?v=dQw4w9WgXcQ", '<p><iframe width="420" height="315" src="//www.youtube.com/embed/dQw4w9WgXcQ" frameborder="0" allowfullscreen></iframe></p>'),

    # Gist link
    ("https://gist.github.com/afrendeiro/6732a46b949e864d6803", '<p><script src="https://gist.github.com/afrendeiro/6732a46b949e864d6803.js"></script></p>'),

    # No linking in code block test
    ("""```print http://www.psu.edu```""", """<p><code>print http://www.psu.edu</code></p>"""),

    # <code> block test
    ("""```print 123```""", """<p><code>print 123</code></p>"""),

    # <pre> <code> block test
    ("""
        print 123
    """, "<pre><code>    print 123</code></pre>"),

    # Test <, >, and & rendering
    ("1 > 0    1 < 2    foo & bar", "<p>1 &gt; 0    1 &lt; 2    foo &amp; bar</p>"),

    # Test <, >, and & in code block
    ("```1 > 0    1 < 2    foo & bar```", "<p><code>1 &gt; 0    1 &lt; 2    foo &amp; bar</code></p>"),

    # Mentioned user
    ("@test", '<p><a href="/accounts/profile/5/" rel="nofollow">tested2</a></p>'),

    # Url auto-link tests
    ("http://www.psu.edu", '<p><a href="http://www.psu.edu" rel="nofollow">http://www.psu.edu</a></p>'),

    # Auto-link enclosed by parenthesis
    ("(http://www.psu.edu)", '<p>(<a href="http://www.psu.edu" rel="nofollow">http://www.psu.edu</a>)</p>'),

    # Test unclosed tags
    ("<b> foo", '<p><b> foo</b></p><b></b>'),

    # Test nested unclosed tags
    ("<b><b><b><b>foo  ", '<p><b><b><b><b>foo</b></b></b></b></p><b><b><b></b></b></b>'),

    # Test ftp
    ("Here is an [ftp link](ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz).", '<p>Here is an <a href="ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz" rel="nofollow">ftp link</a>.</p>')

]


class MarkdownTest(TestCase):
    def setUp(self):
        # Create user
        logger.setLevel(logging.WARNING)
        self.owner = User.objects.create(username="test", email="tested2@tested.com", password="tested")

        self.owner.profile.uid = "5"
        self.owner.profile.handle = "test"
        self.owner.profile.save()
        self.owner.save()

        # Create posts used for anchors
        self.post = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                               type=models.Post.QUESTION, uid="1")

        self.answer = models.Post.objects.create(title="Test", author=self.owner, content="Test",
                                                 type=models.Post.ANSWER, uid="2")

        pass

    def test_markdown(self):

        error_count = 0
        for test in TEST_CASES:
            given, expected = test

            html = markdown.parse(given, clean=True, escape=False)
            html = html.replace("\n", "")

            if expected != html:
                print("\n\n")
                print("--- Input markdown ---")
                print (given)
                print("--- Expected html ---")
                print(expected)
                print("--- Observed html ---")
                print(html)
                error_count += 1

        # Catch all errors at once.
        self.assertTrue(error_count == 0)
