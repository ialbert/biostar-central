from django.test import TestCase
from biostar3.forum.mailer import EmailTemplate
from django.core import mail
import logging

logging.disable(logging.ERROR)

# Create your tests here.
class EmailTests(TestCase):
    def test_emailtemplate(self):
        "Test processing email templates"
        EQ = self.assertEqual
        data = dict(name="John Doe", content="world")
        em = EmailTemplate("mailer_test.html", data=data)

        EQ(em.subj, "Hello John Doe!")
        EQ(em.text, "Hello world!")
        EQ(em.html, "Hello <b>world!</b>")

        em.send(from_email="foo@com", to=["bar@com"])

        # Test that the message has been sent
        EQ(len(mail.outbox), 1)
        EQ(mail.outbox[0].subject, em.subj)
