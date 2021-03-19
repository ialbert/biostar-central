from allauth.socialaccount.adapter import DefaultSocialAccountAdapter
from .models import User, Profile


class SocialAccountAdapter(DefaultSocialAccountAdapter):

    def pre_social_login(self, request, sociallogin):
        """
        Check to see is an account with given email exists in biostars, then logs in.

        User already authenticated with social login at this point.

        """
        if sociallogin.is_existing:
            return

        user = sociallogin.user
        user = User.objects.filter(email=user.email).first()

        if user:
            sociallogin.connect(request, user)
            Profile.objects.filter(user=user).update(email_verified=True)

    def save_user(self, request, sociallogin, form=None):
        """
        Instantiates a new User instance.
        """
        user = super(SocialAccountAdapter, self).save_user(request, sociallogin, form)
        user = User.objects.filter(email=user.email).first()

        # Verify the user email once
        Profile.objects.filter(user=user).update(email_verified=True)


