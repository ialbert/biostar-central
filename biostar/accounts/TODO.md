

Occasional bug, but it happens, invalid password set, fix:


https://stackoverflow.com/questions/62818777/force-django-to-send-reset-password-email-even-if-usable-password-is-not-set-for

def get_users(self, email):
    active_users = get_user_model()._default_manager.filter(
        email__iexact=email, is_active=True)
    return active_users
