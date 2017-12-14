import os

from django.test import  RequestFactory
from django.contrib.messages.storage import fallback

from biostar.engine import tasks


def fake_request(url, data, user, method="POST"):
    "Make a fake request; defaults to POST."

    methods = {"POST": RequestFactory().post, "GET": RequestFactory().get}

    assert method in methods

    request = methods[method](url, data)

    # Mimic messaging system
    request.session = {}
    messages = fallback.FallbackStorage(request=request)
    request._messages = messages

    request.user = user

    return request


def remove_test_folders(target):
    "Remove project/job folder created in media dir when testing "

    target = os.path.abspath(target)

    #tasks.execute(f"rm -rf {target}")
    pass