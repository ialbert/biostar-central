from django.test import TestCase
from biostar3.forum.mailer import EmailTemplate
from django.core import mail

# Create your tests here.
class EmailTests(TestCase):
    def test_emailtemplate(self):
        "Test processing email templates"
        EQ = self.assertEqual
        data = dict(name="World", content="content")
        em = EmailTemplate("mailer_test.html", data=data)

        EQ(em.subj, "Hello World!")
        EQ(em.text, "Text content")
        EQ(em.html, "Html content")

        em.send(from_email="foo@com", to=["bar@com"])

        # Test that the message has been sent
        EQ(len(mail.outbox), 1)
        EQ(mail.outbox[0].subject, em.subj)
