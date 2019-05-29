import hjson
import bleach
from urllib.request import urlopen, Request
import logging
from django.utils.encoding import force_bytes
from django.utils.http import urlsafe_base64_encode
from django.contrib import auth
from django.conf import settings
from django.template import loader

from biostar.emailer.auth import notify
from .models import User, Profile, Message
from . import util, const
from .tokens import account_verification_token

logger = logging.getLogger('engine')


def create_local_messages(template, sender, rec_list, context={}, subject="", uid=None):
    """
    Create batch message from sender for a given recipient_list
    """
    tmpl = loader.get_template(template_name=template)
    html = tmpl.render(context)
    body = bleach.clean(html)

    msgs = []
    for rec in rec_list:
        actual_uid = uid or util.get_uuid(10)
        sent_date = util.now()
        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject, sent_date=sent_date,
                                     uid=actual_uid, body=body, html=html)

        msgs.append(msg)

    return msgs


def check_user(email, password):
    """
    Used to validate user across apps. Returns a tuple ( login message, False or True )
    """

    user = User.objects.filter(email__iexact=email).order_by('-id').first()

    if not user:
        return "This email does not exist.", False

    user = auth.authenticate(username=user.username, password=password)

    if not user:
        return "Invalid Password", False

    if user.profile.state in [Profile.BANNED, Profile.SUSPENDED]:
        msg = f"Login not allowed. Account is <b> {user.profile.get_state_display()}</b>"
        return msg,  False

    elif user and not user.is_active:
        return "This user is not active.", False

    elif user and user.is_active:
        return "Login successful!", True

    return "Invalid fallthrough", False




def send_verification_email(user):

    from_email = settings.DEFAULT_FROM_EMAIL
    userid = urlsafe_base64_encode(force_bytes(user.pk))
    token = account_verification_token.make_token(user)
    template = "accounts/email_verify.html"
    email_list = [user.email]
    context = dict(token=token, userid=userid, user=user)

    # Send the verification email
    notify(template_name=template, email_list=email_list,
           extra_context=context, from_email=from_email,
           subject="Verify your email", send=True)

    return True
