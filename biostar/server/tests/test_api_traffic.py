from datetime import datetime
import json

from django.test import SimpleTestCase
from django.core.urlresolvers import reverse

from ..api import datetime_to_iso, datetime_to_unix


class ApiTrafficTest(SimpleTestCase):
    def test_no_traffic(self):
        """
        There is no posts in the db.
        """
        r = self.client.get(reverse('api-traffic'))
        now = datetime.now()
        content = json.loads(r.content)
        self.assertTrue(content['date'].startswith(datetime_to_iso(now)[:11]))
        self.assertTrue(datetime_to_unix(now) - content['timestamp'] < 100)
        self.assertTrue(content['post_views_last_60_min'] == 0)