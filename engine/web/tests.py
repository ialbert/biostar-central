from django.contrib.auth.models import AnonymousUser, User
from django.test import TestCase, RequestFactory

from .views import custom_login, signup

'''
 if not user:
        # This user needs to be created/
        created = True
        user = User.objects.create(email=email, username=get_uuid())
        user.set_password(get_uuid())
        user.save()
        user.profile.name = name
'''


def get_uuid(limit=32):
    return str(uuid.uuid4())[:limit]


test_urls = [] 

class AuthTestCase(TestCase):

 
    def xtest_signup(self):

        c = Client()
        response = c.post('/signup/', {'email': 'test2@test.come', 'password1': 'test',
                   'password2':'test'})
        print(response.status_code, 'signup')

    def xtest_login(self):

        c = Client()
        response = c.post('/login/', {'username': 'test2@test.come', 'password': 'test'})
        print(response.status_code, 'login')
    




