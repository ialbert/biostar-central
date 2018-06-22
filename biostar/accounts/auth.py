import hjson
import logging

from django.contrib import auth
from .models import User, Profile
from . import util

logger = logging.getLogger('engine')


def check_user(email, password):
    "Used to validate user across apps. Returns a tuple ( login message, False or True ) "

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


def create_user_from_json(json_dict):

    email = json_dict.get("email")
    password = json_dict.get("password", str(util.get_uuid(16)))

    if not email:
        return

    username = util.get_uuid(16)
    name = json_dict.get("name")
    user = User.objects.filter(email=email)

    # Create a user if email is unique
    if not user:
        user = User.objects.create(email=email, username=username, first_name=name,
                                   password=password)
    else:
        logger.error("User with same email already exists.")
        return user.first()

    # Recreate profile
    state = json_dict.get("status", Profile.NEW)
    role = json_dict.get("type", Profile.NORMAL)

    # Given user-id is going to be loaded as a uid
    uid = json_dict.get("id", util.get_uuid(16))
    Profile.objects.filter(user=user).update(uid=uid, state=state, role=role,
                                             website=json_dict.get("website", ""),
                                             scholar=json_dict.get("scholar", ""),
                                             text=json_dict.get("text", ""),
                                             twitter=json_dict.get("twitter", ""),
                                             score=json_dict.get("score", 0),
                                             last_login=json_dict.get("last_login"),
                                             location=json_dict.get("location", ""),
                                             date_joined=json_dict.get("date_joined")
                                             )
    return user