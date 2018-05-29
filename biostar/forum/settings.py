from biostar.settings import *


ALLOWED_TAGS = "p div br code pre h1 h2 h3 h4 hr span s sub sup b i img strong strike em underline super table "
ALLOWED_TAGS += "thead tr th td tbody"
ALLOWED_TAGS = ALLOWED_TAGS.split()

ALLOWED_STYLES = 'color font-weight background-color width height'.split()
ALLOWED_ATTRIBUTES = {
    '*': ['class', 'style'],
    'a': ['href', 'rel'],
    'img': ['src', 'alt', 'width', 'height'],
    'table': ['border', 'cellpadding', 'cellspacing'],

}


# Configure language detection
LANGUAGE_DETECTION = ['en']
