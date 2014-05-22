from django.test import SimpleTestCase
from django.core.urlresolvers import reverse


class ApiTest1(SimpleTestCase):
    def test_no_posts(self):
        """
        There is no posts in the db.
        """
        r = self.client.get(reverse('api-stats-on-day', kwargs={'day': 0}))
        self.assertEqual(r.content, '{}')