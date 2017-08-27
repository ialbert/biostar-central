from django.contrib.auth.models import BaseUserManager
from django.utils import timezone
from django.utils.translation import ugettext_lazy as helpers

class UserManager(BaseUserManager):

    def _create_user(self, email, password, **extra_fields):

        now = timezone.now()

        if not email:
            raise ValueError(helpers(u'The given username must be set'))
        email = self.normalize_email(email)
        user = self.model(email=email,
                          is_active=True,
                          last_login=now,
                          date_joined=now, **extra_fields)
        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_user(self, email, password=None, **extra_fields):
        return self._create_user(email, password, **extra_fields)

