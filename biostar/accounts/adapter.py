from allauth.socialaccount.adapter import DefaultSocialAccountAdapter
from .models import User


class SocialAccountAdapter(DefaultSocialAccountAdapter):

    def pre_social_login(self, request, sociallogin):
        """
        Check to see is an account with given email exists in biostars, then logs in.

        User already authenticated with social login at this point.

        """
        if sociallogin.is_existing:
            return

        user = sociallogin.user
        user = User.objects.filter(email=user.email)

        if user.exists():
            sociallogin.connect(request, user.first())
