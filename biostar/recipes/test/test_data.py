import logging
import os
from unittest.mock import patch, MagicMock

from django.conf import settings
from django.test import TestCase, override_settings
from django.urls import reverse

from biostar.recipes import models, views, auth, const
from biostar.utils.helpers import fake_request, get_uuid

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'tested'))

__CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

logger = logging.getLogger('engine')


@override_settings(MEDIA_ROOT=TEST_ROOT)
class DataViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username=f"tested{get_uuid(10)}", email="tested@l.com")
        self.owner.set_password("tested")

        self.project = auth.create_project(user=self.owner, name="tested", text="Text", summary="summary",
                                           uid="tested")
        self.project.save()

        # Set up generic data for editing
        self.data = auth.create_data(project=self.project, path=__file__, name="tested")

    @patch('biostar.recipes.models.Data.save', MagicMock(name="save"))
    def test_data_edit(self):
        "Test Data edit view with POST request"

        data = {'name': "new_data", 'summary': "summary", 'text': "tested", }

        url = reverse('data_edit', kwargs=dict(uid=self.data.uid))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.data_edit(request=request, uid=self.data.uid)

        obj = {}
        self.data.fill_dict(obj=obj)

        self.assertTrue("toc" in obj, "Table of content not added during fill_dict()")

        self.process_response(response=response, data=data, save=True)

    @patch('biostar.recipes.models.Data.save', MagicMock(name="save"))
    def test_data_upload(self):
        "Test Data upload POST request"

        data = {
            'file': open(__file__, 'r'),
            'summary': 'summary',
            "text": "tested",
        }

        url = reverse('data_upload', kwargs=dict(uid=self.project.uid))

        # Create a new user and give them upload access
        user = models.User.objects.create_user(username="test2", email="test2@l.com")
        user.set_password("tested")
        user.save()
        access = models.Access(access=models.Access.WRITE_ACCESS,
                               user=user,
                               project=self.project)
        access.save()

        request = fake_request(url=url, data=data, user=user)
        response = views.data_upload(request=request, uid=self.project.uid)

        self.process_response(response=response, data=data, save=True)

    def test_add_data(self):
        "Test adding data directory to a project using management commands "

        data_directory = auth.join(__file__, "..", "data")

        data = auth.create_data(project=self.project, path=data_directory)

        self.assertTrue(os.path.exists(data.get_data_dir()), "Directory not being linked")

    def test_data_copy_paste(self):
        "Test data copy and paste interface"

        url = reverse('data_copy', kwargs=dict(uid=self.data.uid))
        clear_url = reverse('clear_clipboard', kwargs=dict(uid=self.project.uid))
        paste_url = reverse('data_paste', kwargs=dict(uid=self.project.uid))
        data = {settings.CLIPBOARD_NAME: const.COPIED_DATA}

        request = fake_request(url=url, data={}, user=self.owner)
        response = views.data_copy(request=request, uid=self.data.uid)
        clear_request = fake_request(url=clear_url, data=data, user=self.owner)
        clear_response = views.clear_clipboard(request=clear_request, uid=self.project.uid)

        self.process_response(response=response, data={})
        self.process_response(response=clear_response, data={})

        # Copy again and paste this time
        copy_request = fake_request(url=url, data={}, user=self.owner)
        copy_response = views.data_copy(request=copy_request, uid=self.data.uid)
        paste_request = fake_request(url=paste_url, data={}, user=self.owner)
        paste_response = views.data_paste(request=paste_request, uid=self.project.uid)

        self.process_response(response=copy_response, data={})
        self.process_response(response=paste_response, data={})

    def test_data_delete(self):
        "Test the data delete"

        url = reverse('data_delete', kwargs=dict(uid=self.data.uid))

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.data_delete(request=request, uid=self.data.uid)

        self.process_response(response=response, data={})

    def test_data_file_copy_paste(self):
        "Test the data file copying interface"
        url = reverse("data_file_copy", kwargs=dict(uid=self.data.uid, path=self.data.get_data_dir()))
        data = {settings.CLIPBOARD_NAME: const.COPIED_FILES, 'path':self.data.get_data_dir()}
        request = fake_request(url=url, data=data, user=self.owner)

        response = views.data_file_copy(request=request, uid=self.data.uid, path=self.data.get_data_dir())

        self.process_response(response, data={})

    def test_data_serve(self):
        "Test data file serving"

        from django.http.response import FileResponse

        data = {"paths": self.data.get_files()[0]}
        url = reverse('data_serve', kwargs=dict(uid=self.data.uid, path=data["paths"]))

        request = fake_request(url=url, data=data, user=self.owner)

        response = views.data_serve(request=request, uid=self.data.uid, path=data["paths"])

        self.assertTrue(isinstance(response, FileResponse), "Response is not a FileResponse type.")

    def process_response(self, response, data, save=False):
        "Check the response on POST request is redirected"

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to project view after tested :\nresponse:{response}")

        if save:
            self.assertTrue(models.Data.save.called, "save() method not called")
