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


ALL_PAGES = ["/", "/projects", "/login", "/signup"]


class SimplePageResponses(TestCase):

    # Test given pages ( with new client for every page )
    def test_pages(self):

        for page in ALL_PAGES:
            client = Client()
            response = client.get(page)
            try:
                self.assertEqual(response.status_code, 200)
            except AssertionError:
                # signup has a redirect response
                self.assertEqual(response.status_code, 301)



    




