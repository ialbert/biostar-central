import logging
from django.contrib.sites.models import Site
from biostar.emailer import sender
from django.core import management
from django.test import TestCase, RequestFactory
from django.conf import settings
from mailer.engine import send_all

logger = logging.getLogger('engine')


class SendMailTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

    def test_send_mail(self):
        "Test email sending using management command."
        management.call_command('test_email')


class ModelTests(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)


    def test_email_group(self):
        from biostar.emailer.models import EmailGroup

        test = EmailGroup(name="test", uid="test")
        test.save()

        print(test)

    def test_email_address(self):
        from biostar.emailer.models import EmailAddress

        test = EmailAddress(name="test", uid="test")
        test.save()

        print(test)
        pass


