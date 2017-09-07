from django.test import TestCase
from django.test import Client

'''
 if not user:
        # This user needs to be created/
        created = True
        user = User.objects.create(email=email, username=get_uuid())
        user.set_password(get_uuid())
        user.save()
        user.profile.name = name
'''



class AuthTestCase(TestCase):


    def test_signup(self):

        c = Client()
        response = c.post('/signup/', {'email': 'test2@test.come', 'password1': 'test',
                   'password2':'test'})
        print(response.status_code, 'signup')

    def test_login(self):

        c = Client()
        response = c.post('/login/', {'username': 'test@test.come', 'password': 'test'})
        print(response.status_code, 'login')




