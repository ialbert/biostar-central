import logging
from django.core import management
from django.test import TestCase

from django.core.files import File
from django.urls import reverse
from biostar.engine import models, views, auth, factory, forms
from biostar.engine import util as engine_util

from . import util

logger = logging.getLogger('engine')


class SiteAdminTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.user = models.User.objects.create_superuser(username="test", email="test@test.com",
                                                         password="test")
        self.user.save()

    def test_site_admin(self):
        "Test site admin page"

        url =  reverse("site_admin")
        request = util.fake_request(url=url, data={}, method="GET",user=self.user)

        response = views.site_admin(request=request)
        self.assertEqual(response.status_code, 200, "Can not load site admin page specific to biostar-engine.")


class FactoryTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)


        owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=owner, name="test",
                                           text="Text", summary="summary", uid="testing")
        self.json_data = {
                "label":"test",
                "help":"Test json data",
                "choices":["test1", "test2"]
        }

    def test_factory_fields(self):
        "Testing factory module that generates fields"

        util.remove_test_folders(self.project.get_project_dir())

        # All valid field types.
        field_types = factory.get_field_types()

        for display_type in field_types:

            # Test that each field type can be rendered.
            json_data = dict(display=display_type)

            field = factory.dynamic_field(json_data)
            if not field:
                message = f"field generator for display={display_type} failed"
                self.assertFalse(message)


    def test_dynamic_field(self):
        "Test data generator"

        from biostar.engine import const

        pre = models.Data.objects.count()
        data = auth.create_data(self.project, path=__file__)
        post = models.Data.objects.count()

        util.remove_test_folders(self.project.get_project_dir())

        self.assertTrue(post == (pre + 1), "Error creating data in database")

        display_type = const.DROPDOWN

        json_data = dict(display=display_type, path=data.get_path())

        field = factory.dynamic_field(json_data, project=self.project)

        if not field:
            message = f"field generator for display={display_type} failed"
            self.assertFalse(message)


    def test_data_generator(self):
        "Test data generator"
        util.remove_test_folders(self.project.get_project_dir())
        pass


class ManagementCommandTest(TestCase):

    def setUp(self):

        logger.setLevel(logging.WARNING)

        self.owner = models.User.objects.filter(is_superuser=True).first()

        self.project = auth.create_project(user=self.owner, name="test",
                                           text="Text", summary="summary", uid="testing")
        self.analysis = auth.create_analysis(project=self.project, json_text='{test:{value:"test"}}',
                                             template="echo {{test.value}}", security=models.Analysis.AUTHORIZED)
        self.analysis.save()
        self.job = auth.create_job(analysis=self.analysis)
        self.job.save()


    def test_add_data(self):
        "Test adding data to a project using management commands "

        pre = models.Data.objects.all().count()
        management.call_command('data', path=__file__, uid="testing")
        post = models.Data.objects.all().count()

        util.remove_test_folders(self.project.get_project_dir())
        util.remove_test_folders(self.job.path)

        self.assertTrue(post == (pre + 1), "Error creating adding in database with management command")


    def test_job_runner(self):
        "Testing Job runner using management command"

        management.call_command('job', id=self.job.id, verbosity=2, list=True)

        util.remove_test_folders(self.project.get_project_dir())
        util.remove_test_folders(self.job.path)

    def test_job_commands(self):
        "Testing extra job runner commands"

        management.call_command('job', verbosity=2, list=True)


class UtilTests(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

    def test_smart_preview(self):

        collect = auth.findfiles("biostar/engine/test/data", collect=[])
        for fname in collect:
            text = engine_util.smart_preview(fname)
            self.assertTrue("error" not in text.split(), f"Preview error with {fname}")
