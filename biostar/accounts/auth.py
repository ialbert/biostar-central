import hjson

from django.contrib import auth
from .models import User, Profile


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
        return  "Login successful!", True


    return "Invalid fallthrough", False



def create_user_from_json(json_dict):

    user = ''


    print(json_dict)
    1/0




    return user