# django-openid-auth -  OpenID integration for django.contrib.auth
#
# Copyright (C) 2010 Canonical Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import unittest

from django.test import TestCase

from django_openid_auth.auth import OpenIDBackend
from openid.consumer.consumer import SuccessResponse
from openid.consumer.discover import OpenIDServiceEndpoint
from openid.message import Message, OPENID2_NS


SREG_NS = "http://openid.net/sreg/1.0"
AX_NS = "http://openid.net/srv/ax/1.0"

class OpenIDBackendTests(TestCase):

    def setUp(self):
        super(OpenIDBackendTests, self).setUp()
        self.backend = OpenIDBackend()

    def test_extract_user_details_sreg(self):
        endpoint = OpenIDServiceEndpoint()
        message = Message(OPENID2_NS)
        message.setArg(SREG_NS, "nickname", "someuser")
        message.setArg(SREG_NS, "fullname", "Some User")
        message.setArg(SREG_NS, "email", "foo@example.com")
        response = SuccessResponse(
            endpoint, message, signed_fields=message.toPostArgs().keys())

        data = self.backend._extract_user_details(response)
        self.assertEqual(data, {"nickname": "someuser",
                                "first_name": "Some",
                                "last_name": "User",
                                "email": "foo@example.com"})

    def test_extract_user_details_ax(self):
        endpoint = OpenIDServiceEndpoint()
        message = Message(OPENID2_NS)
        attributes = [
            ("nickname", "http://axschema.org/namePerson/friendly", "someuser"),
            ("fullname", "http://axschema.org/namePerson", "Some User"),
            ("email", "http://axschema.org/contact/email", "foo@example.com"),
            ]
        message.setArg(AX_NS, "mode", "fetch_response")
        for (alias, uri, value) in attributes:
            message.setArg(AX_NS, "type.%s" % alias, uri)
            message.setArg(AX_NS, "value.%s" % alias, value)
        response = SuccessResponse(
            endpoint, message, signed_fields=message.toPostArgs().keys())

        data = self.backend._extract_user_details(response)
        self.assertEqual(data, {"nickname": "someuser",
                                "first_name": "Some",
                                "last_name": "User",
                                "email": "foo@example.com"})

    def test_extract_user_details_ax_split_name(self):
        endpoint = OpenIDServiceEndpoint()
        message = Message(OPENID2_NS)
        attributes = [
            ("nickname", "http://axschema.org/namePerson/friendly", "someuser"),
            # Include this key too to show that the split data takes
            # precedence.
            ("fullname", "http://axschema.org/namePerson", "Bad Data"),
            ("first", "http://axschema.org/namePerson/first", "Some"),
            ("last", "http://axschema.org/namePerson/last", "User"),
            ("email", "http://axschema.org/contact/email", "foo@example.com"),
            ]
        message.setArg(AX_NS, "mode", "fetch_response")
        for (alias, uri, value) in attributes:
            message.setArg(AX_NS, "type.%s" % alias, uri)
            message.setArg(AX_NS, "value.%s" % alias, value)
        response = SuccessResponse(
            endpoint, message, signed_fields=message.toPostArgs().keys())

        data = self.backend._extract_user_details(response)
        self.assertEqual(data, {"nickname": "someuser",
                                "first_name": "Some",
                                "last_name": "User",
                                "email": "foo@example.com"})

    def test_extract_user_details_ax_broken_myopenid(self):
        endpoint = OpenIDServiceEndpoint()
        message = Message(OPENID2_NS)
        attributes = [
            ("nickname", "http://schema.openid.net/namePerson/friendly",
             "someuser"),
            ("fullname", "http://schema.openid.net/namePerson", "Some User"),
            ("email", "http://schema.openid.net/contact/email",
             "foo@example.com"),
            ]
        message.setArg(AX_NS, "mode", "fetch_response")
        for (alias, uri, value) in attributes:
            message.setArg(AX_NS, "type.%s" % alias, uri)
            message.setArg(AX_NS, "value.%s" % alias, value)
        response = SuccessResponse(
            endpoint, message, signed_fields=message.toPostArgs().keys())

        data = self.backend._extract_user_details(response)
        self.assertEqual(data, {"nickname": "someuser",
                                "first_name": "Some",
                                "last_name": "User",
                                "email": "foo@example.com"})

def suite():
    return unittest.TestLoader().loadTestsFromName(__name__)
