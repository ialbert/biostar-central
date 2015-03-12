#
# Celery config file.
#
# Must be imported into the settings


# This must be set to true to route tasks via celery.
CELERY_ENABLED = False

# The serializer format.
CELERY_TASK_SERIALIZER = 'json'

# Ignore other content
CELERY_ACCEPT_CONTENT = ['json']


CELERY_RESULT_SERIALIZER = 'json'

CELERY_ENABLE_UTC = True