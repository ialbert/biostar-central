from django.db import models
from django.contrib.auth.models import AbstractBaseUser


class User(AbstractBaseUser):
    # Cutsom user class based on email 

    email = models.EmailField(verbose_name='email',
        max_length=254, unique=True, db_index=True,
        help_text='Required.'
    )
    
    USERNAME_FIELD = 'email'

    def __str__(self):
        return self.email
        


