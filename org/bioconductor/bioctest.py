from biostar.settings.base import *
import os

# Add the custom directory

# This should be a fully qualified absolute path
THEME_PATH = os.path.join(HOME_DIR, "org", "bioconductor", "templates")

TEMPLATE_DIRS = [THEME_PATH] + TEMPLATE_DIRS
