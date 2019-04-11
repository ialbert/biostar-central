import logging, time, shutil, subprocess
from django.core import management
from django.utils.encoding import force_text


import time

logger = logging.getLogger("engine")

HAS_UWSGI = False


COUNTER = 1

try:
    from uwsgidecorators import *

    HAS_UWSGI = True


    def send_mail():
        return

    def create_user_awards(user):
        return


    def check_user_profile(ip, user):
        return

except ModuleNotFoundError as exc:
    pass