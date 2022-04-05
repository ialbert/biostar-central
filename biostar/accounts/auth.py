import logging
from django.utils.encoding import force_bytes
from django.utils.http import urlsafe_base64_encode
from django.contrib import auth
from django.template import loader
from django.conf import settings


from biostar.emailer.tasks import send_email
from .models import User, Profile
from .tokens import account_verification_token

logger = logging.getLogger('engine')


def validate_login(email, password):
    """
    Used to validate user across apps. Returns a tuple ( login message, False or True )
    """

    user = User.objects.filter(email__iexact=email).order_by('-id').first()

    if not user:
        return "This email does not exist.", False

    user = auth.authenticate(username=user.username, password=password)

    if not user:
        return "Invalid Password", False

    if not user.profile.is_valid:
        msg = f"Login not allowed. Account is not valid."
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
    subject = "accounts/email_verify_subject.md"
    subject = loader.get_template(template_name=subject).render(dict())

    # Send the verification email
    send_email(template_name=template, recipient_list=email_list,
               extra_context=context, from_email=from_email, subject=subject)

    return True

