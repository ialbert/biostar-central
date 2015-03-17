"""
Award specifications
"""
import logging
from .models import Badge, Award

logger = logging.getLogger('biostar')

class AwardDef(object):
    def __init__(self, name, desc, func, icon, type=Badge.BRONZE):
        self.name = name
        self.desc = desc
        self.fun = func
        self.icon = icon
        self.template = "badge/default.html"
        self.type = type

    def validate(self, *args, **kwargs):
        try:
            value = self.fun(*args, **kwargs)
            return value
        except Exception, exc:
            logger.error("validator error %s" % exc)
        return 0

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)