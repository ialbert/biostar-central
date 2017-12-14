from django.test import  RequestFactory
from django.contrib.messages.storage import fallback


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
