import logging
import os
from unittest.mock import patch, MagicMock

from django.test import TestCase, override_settings
from django.urls import reverse

from biostar.engine import models, views, auth
from django.conf import settings


logger = logging.getLogger('engine')


class PostTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        pass








