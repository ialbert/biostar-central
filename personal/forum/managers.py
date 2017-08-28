
from django.contrib.auth.base_user import BaseUserManager

class UserManager(BaseUserManager):

    use_in_migrations = True

    def _create_user(self, email, password, username=None, **extra_fields):
        """
        Create and save a user with the given emails and password.
        """
        if not email:
            raise ValueError('The given email must be set')
        email = self.normalize_email(email)
        #username = self.model.normalize_username(username)
        user = self.model(email=email, **extra_fields)
        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_user(self, email, password, username=None, **extra_fields):
    
        return self._create_user(email, password, 
                                is_staff = False, 
                                is_active = True,
                                 **extra_fields)

    def create_superuser(self, email, password, username=None, **extra_fields):

        return self._create_user(email, password, 
                                is_staff = True, 
                                is_active = True,
                                 **extra_fields)


