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


def create_messages(template, sender, rec_list, context={}, subject=""):
    """
    Create batch message from sender for a given recipient_list
    """
    tmpl = loader.get_template(template_name=template)
    html = tmpl.render(context)
    body = bleach.clean(html)

    msgs = []
    for rec in rec_list:
        msg = Message.objects.create(sender=sender, recipient=rec, subject=subject, body=body, html=html)
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


def check_user_profile(request, user):

    # Get the ip information
    ip1 = request.META.get('REMOTE_ADDR', '')
    ip2 = request.META.get('HTTP_X_FORWARDED_FOR', '').split(",")[0].strip()
    ip = ip1 or ip2 or '0.0.0.0'

    logger.debug(f"profile check from {ip} on {user}")
    # Check and log location.
    if not user.profile.location:
        try:
            url = f"http://api.hostip.info/get_json.php?ip={ip}"
            req = Request(url=url, headers={'User-Agent': 'Mozilla/5.0'})
            user_info = urlopen(req, timeout=3).read()
            logger.debug(f"{ip}, {user}, {url}")

            data = hjson.loads(user_info)
            location = data.get('country_name', '').title()
            location = "localhost" if ip in ('127.0.0.1') else location

            if "unknown" not in location.lower():
                Profile.objects.filter(user=user).update(location=location)
        except Exception as exc:
            logger.error(exc)


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