from django.db import models
from django.contrib.auth.models import AbstractBaseUser
from django.utils.translation import gettext, gettext_lazy as helpers
from .managers import UserManager

class User(AbstractBaseUser):
    # Cutsom user class based on email 

    email = models.EmailField(verbose_name='email',
        max_length=254, unique=True, db_index=True,
        help_text='Required.'
    )

    is_staff = models.BooleanField(
        helpers('staff status'),
        default=False,
        help_text=helpers('Designates whether the user can log into this admin site.'),
    )

    is_active = models.BooleanField(
        helpers('active'),
        default=True,
        help_text=helpers(
            'Designates whether this user should be treated as active. '
            'Unselect this instead of deleting accounts.'
        ),
    )

    objects = UserManager()

    USERNAME_FIELD = 'email'

    def __str__(self):
        return self.email
        


