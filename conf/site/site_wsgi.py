import os, logging
from django.core.wsgi import get_wsgi_application
from django.core.management import call_command

logger= logging.getLogger("biostar")

# Override the DJANGO SETTINGS MODULE

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "conf.site.site_settings")

dj = os.environ["DJANGO_SETTINGS_MODULE"]

print(f"*** DJANGO_SETTINGS_MODULE={dj}")
print('**** Mounting app ')
application = get_wsgi_application()

# Initialize recurring tasks.
try:
    import uwsgidecorators

    @uwsgidecorators.timer(30)
    def send_queued_mail(num):
        """Send queued mail every 30 seconds"""
        pass

    @uwsgidecorators.timer(3600)
    def savedata(num):
        """Save the data every hour"""
        pass

except Exception as exc:
    logger.error(f"uwsgidecorators not enabled, {exc}")
