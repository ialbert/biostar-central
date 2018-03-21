import logging
from django.core import management
from biostar.emailer import auth
from django.test import TestCase
from biostar.emailer import models

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


    def test_add_subs(self):
        "Test adding subscription using auth"

        group = models.EmailGroup()
        group.save()
        auth.add_subscription(email="Test@test.com", group=group, name='testing')


    def test_mailing(self):
        "Test sending email using manage commands."
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


