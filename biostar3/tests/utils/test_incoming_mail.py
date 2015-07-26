from django.test import TestCase

from biostar3.utils.incoming_mail import IncomingMail, IncomingMailException
from biostar3.forum.models import User, ReplyToken
from biostar3.forum import auth

from faker import Factory
f = Factory.create()

class IncomingMailTests(TestCase):
    def test_empty_mail(self):
        with self.assertRaises(IncomingMailException):
            IncomingMail.from_raw('')

    def test_valid_address(self):
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = 'action+token+group@biostar.test'
        )

        IncomingMail(data)

    def test_short_address(self):
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = 'action+token@biostar.test'
        )

        with self.assertRaisesRegexp(IncomingMailException, '\[to_address\]'):
            IncomingMail(data)

    def test_long_address(self):
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = 'action+token+group+extra@biostar.test'
        )

        with self.assertRaisesRegexp(IncomingMailException, '\[to_address\]'):
            IncomingMail(data)

    def test_remove_quoted_with_quote(self):
        reply_text = f.sentence()
        quote_text = f.sentence()
        text_with_quote = "{}\n\n> {}".format(reply_text, quote_text)

        data = dict(
            subject = f.sentence(),
            text = text_with_quote,
            to_address = 'action+token+group@biostar.test'
        )

        mail = IncomingMail(data)

        self.assertRegexpMatches(mail.text, reply_text)
        self.assertRegexpMatches(mail.text, quote_text)

        mail.remove_quoted_text()

        self.assertRegexpMatches(mail.text, reply_text)
        self.assertNotRegexpMatches(mail.text, quote_text)

    def test_remove_quoted_without_quote(self):
        reply_text = f.sentence()

        data = dict(
            subject = f.sentence(),
            text = reply_text,
            to_address = 'action+token+group@biostar.test'
        )

        mail = IncomingMail(data)

        self.assertRegexpMatches(mail.text, reply_text)

        mail.remove_quoted_text()

        self.assertRegexpMatches(mail.text, reply_text)

    def test_new_action_invalid_token(self):
        invalid_token = 'INVALID'
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = '{}+{}+group@biostar.test'.format(IncomingMail.NEW, invalid_token)
        )

        mail = IncomingMail(data)

        with self.assertRaisesRegexp(IncomingMailException, '\[author\]'):
            mail.perform_action()

    def test_new_action_creates_post(self):
        user = User.objects.create(name=f.name(), email=f.email())
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = '{}+{}+group@biostar.test'.format(IncomingMail.NEW, user.profile.uuid)
        )

        mail = IncomingMail(data)
        post = mail.perform_action()

        self.assertEqual(post.author, user)

    def test_reply_action_invalid_token(self):
        invalid_token = 'INVALID'
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = '{}+{}+group@biostar.test'.format(IncomingMail.REPLY, invalid_token)
        )

        mail = IncomingMail(data)

        with self.assertRaisesRegexp(IncomingMailException, '\[reply_token\]'):
            mail.perform_action()

    def test_reply_action_creates_post(self):
        user = User.objects.create(name=f.name(), email=f.email())
        toplevel_post = auth.create_toplevel_post(data=dict(title=f.sentence(), content=f.sentence()), user=user)
        reply_token = ReplyToken.objects.create(user=user, post=toplevel_post)
        data = dict(
            subject = f.sentence(),
            text = f.sentence(),
            to_address = '{}+{}+group@biostar.test'.format(IncomingMail.REPLY, reply_token.token)
        )

        mail = IncomingMail(data)
        post = mail.perform_action()

        self.assertEqual(post.author, user)
        self.assertEqual(post.parent, toplevel_post)