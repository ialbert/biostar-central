import logging

from biostar.emailer import auth
from django.test import TestCase

logger = logging.getLogger('engine')


class SendMailTest(TestCase):


    def setUp(self):
        logger.setLevel(logging.WARNING)

    def test_send_mail(self):
        "Test email sending using auth."

        context = dict(target_email="2@lvh.me")
        from_mail= "mailer@biostars.org"
        template_name = "test_email.html"
        successful = auth.notify(email_list=["2@lvh.me"], extra_context=context,
                                 template_name=template_name, from_email=from_mail)

        self.assertTrue(successful, "Error sending mail")



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


