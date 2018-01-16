import logging
from django.test import TestCase
from unittest.mock import patch, MagicMock
from django.urls import reverse
from django.core import management

from biostar.tools import const
from biostar.engine import models, views, auth
from . import util


logger = logging.getLogger('engine')


class DataViewTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        # Set up generic owner
        self.owner = models.User.objects.create_user(username="test", email="test@l.com")
        self.owner.set_password("test")

        self.project = auth.create_project(user=self.owner, name="test", text="Text", summary="summary",
                                           uid="testing")
        self.project.save()

        # Set up generic data for editing
        pre = models.Data.objects.count()
        self.data = auth.create_data(project=self.project, path=__file__)
        self.assertTrue(models.Data.objects.count() == (pre + 1), "Error creating Data in database")

    def Xtest_data_copy_view(self):

        "Test Data copy (create a new project and copy) in views with POST request"

        # 0 is the option picked whe creating and copying
        data = {"project":0}
        url = reverse('data_view', kwargs=dict(id=self.data.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.data_view(request=request, id=self.data.id)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect data view after copying Data:\nresponse:{response}")


    @patch('biostar.engine.models.Data.save', MagicMock(name="save"))
    def test_data_edit(self):
        "Test Data edit view with POST request"

        data = {'name':"new_data", 'summary':"summary", 'text':"testing",
                'sticky':True}

        url = reverse('data_edit', kwargs=dict(id=self.data.id))

        request = util.fake_request(url=url, data=data, user=self.owner)

        response = views.data_edit(request=request, id=self.data.id)

        obj ={}
        self.data.fill_dict(obj=obj)

        self.assertTrue("toc" in obj, "Table of content not added during fill_dict()")

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to data view after copying Data:\nresponse:{response}")

        self.assertTrue( models.Data.save.called, "data.save() method not called when editing.")


    @patch('biostar.engine.models.Data.save', MagicMock(name="save"))
    def test_data_upload(self):
        "Test Data upload POST request"

        data = {'file':open(__file__, 'r'), 'summary':'summary', "text":"testing", "sticky":True}
        url = reverse('data_upload', kwargs=dict(uid=self.project.uid))

        # Create a new user and give them upload access
        user = models.User.objects.create_user(username="test2", email="test2@l.com")
        user.set_password("test")
        user.save()

        access = models.Access(access=models.Access.UPLOAD_ACCESS,
                              user=user,
                              project=self.project)
        access.save()

        request = util.fake_request(url=url, data=data, user=user)
        response = views.data_upload(request=request, uid=self.project.uid)

        self.assertEqual(response.status_code, 302,
                         f"Could not redirect to after uploading:\nresponse:{response}")

        self.assertTrue( f"/data/list/{self.project.uid}/" == response.url,
                         f"Could not redirect to data list after uploading:\nresponse:{response}")

        self.assertTrue( models.Data.save.called, "data.save() method not called when uploading.")


    def test_add_data(self):
        "Test adding data to a project using management commands "

        pre = models.Data.objects.all().count()
        management.call_command('data', path=__file__, uid="testing")
        post = models.Data.objects.all().count()

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")

    def test_data_download(self):
        "Testing data download view"

        from django.http.response import FileResponse

        url = reverse('data_download', kwargs=dict(id=self.data.id))
        request = util.fake_request(url=url, data={}, user=self.owner, method="GET")

        response = views.data_download(request=request, id=self.data.id)

        self.assertTrue(isinstance(response, FileResponse),
                        "Error returning valid FileResponse when downloading")

