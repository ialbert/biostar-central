from django.db import models
from django.contrib.auth.models import AbstractBaseUser


class User(AbstractBaseUser):

    email = models.EmailField(verbose_name='email',
        max_length=254, unique=True, db_index=True,
    )
    USERNAME_FIELD = 'email'

    @property
    # too hacky?
    def username(self):
        return getattr(self, self.USERNAME_FIELD)

    @username.setter
    def set_username(self, value):
        self.email = value
