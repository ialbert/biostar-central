import logging, os
from unittest.mock import patch, MagicMock
from django.test import Client
from django.contrib.auth.models import AnonymousUser
from django.test import TestCase, override_settings
from django.shortcuts import reverse
from biostar.utils.helpers import fake_request
from biostar.accounts.views import user_login
from biostar.accounts.middleware import limiter

# Set the rate limit really low for testing
RATE = '1/d'

# IP request to add to request META
REMOTE_ADDR = '127.0.0.1'


# Whitelist the IP triplet
WHITELIST_IP = [REMOTE_ADDR[:-2]]


@override_settings(RATELIMIT_RATE=RATE)
class RateLimiterTest(TestCase):

    def setUp(self):
        # Set the IP address and host name in request META
        rmeta = dict(REMOTE_ADDR=REMOTE_ADDR)

        # View to call
        self.url = reverse('login')

        # Create the fake request
        self.request = fake_request(url=self.url, data={}, rmeta=rmeta, user=AnonymousUser())

    def run_one(self):
        response = limiter(user_login)(request=self.request)
        return response

    def limiting(self, offset=0, msg=''):
        """
        Encapsulate rate limiting process.
        """
        # Go over the threshold rate
        rate = int(RATE[0]) + offset
        msg = msg or 'User was rate limited before reaching limit'
        # Make repeated requests and reach rate limit
        for r in range(rate):
            # Hit the view
            response = self.run_one()

            # Ensure user is not rate limited yet.
            self.assertTrue(len(response.content) > 0, msg)

    def test_banning(self):
        """
        Test limiter is functioning
        """
        # Iterate and call view until
        # limiting threshold is reached.
        self.limiting()

        # Call one more time to ensure that rate limiter is working
        response = self.run_one()
        self.assertTrue(len(response.content) == 0, 'User was not rate limited at all.')

    @override_settings(WHITELIST_IP=WHITELIST_IP)
    def test_whitelisting_ip(self):
        """
        Test IP whitelisting
        """
        # Iterate and call view until
        # num requests pass threshold by offset.
        self.limiting(offset=10, msg='User was rate limited and has a whitelisted ip.')
