from django.core.urlresolvers import reverse, resolve
from django.test import TestCase
from django.test import Client
from django.core import management
from biostar.engine.models import User
from .urls import urlpatterns



def generate_full_url(pattern):
    return




class GeneralTestCase(TestCase):

    def test_page_response(self):

        #response = self.client.get(reverse('create', args=[self.userName]))

        # Get a super user.
        super_user = User.objects.filter(is_superuser=True).first()
        self.assertTrue(isinstance(super_user, User))

        # Make an admin project with one analysis and job.

        project_id = management.call_command('project', next=True, creator_email=super_user.email)

        #self.assertEqual(response.status_code, 200)
        for pattern in urlpatterns:
            print(pattern)

        self.assertEqual(1+1,4)

        return


    pass










