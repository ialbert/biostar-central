import logging
from pyzmail import PyzMessage

from biostar3.utils.email_reply_parser import EmailReplyParser
from biostar3.forum.models import User, ReplyToken
from biostar3.forum import auth

logger = logging.getLogger("biostar")

class IncomingMailException(Exception):
    def __init__(self, kind, msg):
        self.kind = kind
        self.msg = msg
    def __str__(self):
        return "Incoming Mail Error [{}]: {}".format(self.kind, self.msg)

class IncomingMail:
    NEW = "new"
    REPLY = "reply"
    ACTIONS = [NEW, REPLY]

    @staticmethod
    def from_raw(content):
        data = dict()

        try:
            msg = PyzMessage.factory(content)
        except Exception as exc:
            logger.error(exc)
            content = content.encode('utf8', errors='ignore')
            msg = PyzMessage.factory(content)

        try:
            data['subject'] = msg.get_subject()
        except Exception as exc:
            raise IncomingMailException("subject", exc)

        try:
            data['to_address'] = msg.get_addresses('to')[0][1]
        except Exception as exc:
            raise IncomingMailException("to_address", exc)

        try:
            part = msg.text_part or msg.html_part
            enc = part.charset
            data['text'] = part.get_payload().decode(enc)
        except Exception as exc:
            raise IncomingMailException("text", exc)

        return IncomingMail(data)

    def __init__(self, data):
        self.subject = data.get('subject')

        try:
            self.action, self.token, self.group = data.get('to_address').split('@')[0].split('+')
        except Exception as exc:
            raise IncomingMailException("to_address", exc)

        self.text = data.get('text')

    def remove_quoted_text(self):
        self.text = EmailReplyParser.parse_reply(self.text)

    def perform_action(self):
        if not self.action in IncomingMail.ACTIONS:
            raise IncomingMailException("action", self.action)

        if self.action == IncomingMail.NEW:
            # The token must match a users uuid.
            author = User.objects.filter(profile__uuid=self.token).first()
            if not author:
                raise IncomingMailException("author", "{} does not exist".format(self.token))

            # Create the post in the usergroup
            data = dict(title=self.subject, content=self.text, tags="via email")
            return auth.create_toplevel_post(data=data, user=author)
        elif self.action == IncomingMail.REPLY:
            reply_token = ReplyToken.objects.filter(token=self.token).first()
            if not reply_token:
                raise IncomingMailException("reply_token", "{} does not exist".format(self.token))

            return auth.create_content_post(content=self.text, parent=reply_token.post, user=reply_token.user)