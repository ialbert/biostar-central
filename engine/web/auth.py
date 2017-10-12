# Access authorization to different objects

from django.conf import settings
from engine.models import *


def get_data(obj_id):
    # TODO verify that user may access the data
    data = Data.objects.filter(id=obj_id).first()
    return data

