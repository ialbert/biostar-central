import logging
import os
from django.core import management
from django.test import TestCase, override_settings
from django.forms import ValidationError
from django.core.files import File
from django.urls import reverse
from biostar.recipes import models, views, auth, factory, forms, const, api
from biostar.recipes import util as engine_util
from django.conf import settings


from biostar.utils.helpers import fake_request, get_uuid

TEST_ROOT = os.path.abspath(os.path.join(settings.BASE_DIR, 'export', 'test'))
TOC_ROOT = os.path.join(TEST_ROOT, 'toc')
logger = logging.getLogger('engine')

# Ensure that the table of directory exists.
os.makedirs(TOC_ROOT, exist_ok=True)

class Bunch(object):
    last_valid = template = ''

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)



@override_settings(MEDIA_ROOT=TEST_ROOT, TOC_ROOT=TOC_ROOT)
class SiteAdminTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

        self.user = models.User.objects.create_superuser(username=f"tested{get_uuid(10)}", email="tested@tested.com",
                                                         password="tested")
        self.user.save()

    def test_site_admin(self):
        "Test site admin page"

        url = reverse("site_admin")
        request = fake_request(url=url, data={}, method="GET",user=self.user)

        response = views.site_admin(request=request)
        # admin page specific to biostar-engine.
        self.assertEqual(response.status_code, 200, "Can not load admin page")

    def test_bin_view(self):
        "Test recycle bin view"

        url = reverse('recycle_bin')
        request = fake_request(url=url, data={}, method="GET",user=self.user)
        response = views.recycle_bin(request=request)
        self.assertEqual(response.status_code, 200, "Can not load recyle bin")


@override_settings(MEDIA_ROOT=TEST_ROOT , MULTI_THREAD=False)
class FactoryTest(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)
        self.owner = models.User.objects.filter(is_superuser=True).first()
        self.project = auth.create_project(user=self.owner, name="tested",
                                           text="Text", summary="summary", uid="tested")

    def Xtest_api_change_obj(self):
        """
        Change object image
        """

        new_img = open(os.path.join(TEST_ROOT, "data", "image.png"))

        api.change_image(obj=self.project, fname=new_img)


        return

    def test_factory_fields(self):
        "Testing factory module that generates fields"

        # All valid field types.
        field_types = factory.get_field_types()

        for display_type in field_types:

            # Test that each field type can be rendered.
            json_data = dict(display=display_type)

            field = factory.dynamic_field(json_data)
            if not field:
                message = f"field generator for display={display_type} failed"
                if display_type == const.SQL:
                    return
                self.assertFalse(message)

    def test_dynamic_field(self):
        "Test data generator"

        from biostar.recipes import const

        data = auth.create_data(self.project, path=__file__)

        display_type = const.DROPDOWN

        json_data = dict(display=display_type, value=data.name,
                         source= 'PROJECT')

        field = factory.dynamic_field(json_data, project=self.project)

        if not field:
            self.assertFalse(f"field generator for display={display_type} failed")

    def test_data_generator(self):
        "Test data generator"

        field = factory.data_field_generator(field={}, type="DATA", project=self.project)

        if not field:
            self.assertFalse(f"data field generator failed")

    def test_import_file(self):
        "Test import files tab view"
        url = reverse('root_list')

        request = fake_request(url=url, data={}, user=self.owner)

        response = views.import_files(request)

        self.assertEqual(response.status_code, 200, f"Error with file listing in import tab.")


@override_settings(MEDIA_ROOT=TEST_ROOT)
class UtilTests(TestCase):

    def setUp(self):
        logger.setLevel(logging.WARNING)

    def test_smart_preview(self):

        collect = engine_util.findfiles("biostar/recipes/test/data", collect=[])
        for fname in collect:

            text = engine_util.smart_preview(fname)
            self.assertTrue("error" not in text.split(), f"Preview error with {fname}")

            with self.assertRaises(ValidationError):
                forms.check_size(File(open(fname, "r")), maxsize=0.000001)

