import logging
import os
import shutil
from django.core import management
from django.urls import reverse
from django.test import TestCase, override_settings
from django.conf import settings
from biostar.forum import models, views
from biostar.utils.helpers import fake_request
from biostar.accounts.models import User
from biostar.utils import spamlib

logger = logging.getLogger('engine')

TEST_SPAM_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'test', 'test_spammers'))
TEST_SPAM_DIR = TEST_SPAM_ROOT
TEST_SPAM_INDEX_NAME = "test_spam"


@override_settings(SPAM_INDEX_NAME=TEST_SPAM_INDEX_NAME,
                   SPAM_INDEX_DIR=TEST_SPAM_DIR,
                   CLASSIFY_SPAM=True)
class TestSpam(TestCase):

    def setUp(self):

        # Add test post as spam
        self.owner = User.objects.create(username=f"test", email="tested@tested.com", password="tested",
                                         is_superuser=False, is_staff=False)
        for r in range(10):
            content = f"this is spam {r}"
            self.spam = models.Post.objects.create(title=f"Test spam", author=self.owner,
                                                   content=content, spam=models.Post.SPAM,
                                                   type=models.Post.QUESTION)

        #spam.build_spam_index()

    def Xtest_score(self):
        """
        Test spam scoring
        """
        # Add this post to the spam index.
        spam.add_spam(post=self.spam)

        content = "this is spam"
        self.new_spam = models.Post.objects.create(title=f"Test spam", author=self.owner,
                                               content=content, spam=models.Post.SPAM,
                                               type=models.Post.QUESTION)

        spam.score(post=self.new_spam, threshold=0.0)

        new_spam = models.Post.objects.filter(uid=self.new_spam.uid).first()

        self.assertTrue(new_spam.is_spam, "Spam is classifier is not working")

        pass
