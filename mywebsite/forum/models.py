from django.db import models
from django.contrib.auth.models import AbstractBaseUser

class User(AbstractBaseUser):

    email = models.EmailField(verbose_name='email address',
        max_length=254, unique=True, db_index=True,
    )
    USERNAME_FIELD = 'email'


