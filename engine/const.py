import os
from biostar.tools.const import *

__THIS = os.path.dirname(__file__)


def join(*args):
    return os.path.abspath(os.path.join(*args))


HOME_ICON = "home"
PROJECT_LIST_ICON = "database"
PROJECT_ICON = "archive"
DATA_LIST_ICON = "file text"
DATA_ICON = "file"
ANALYSIS_LIST_ICON = "settings"
ANALYSIS_ICON = "setting"
RESULT_LIST_ICON = "bar chart"
RESULT_ICON = "line chart"
RESULT_VIEW_ICON = "line chart"
RESULT_INDEX_ICON = "list layout icon"
LOGIN_ICON = "sign in"
LOGOUT_ICON = "sign out"
INFO_ICON = "info circle icon"
SIGNUP_ICON = "add user icon"
USER_ICON = "user icon"



FIELD_VISIBLE = "visible"
FIELD_ORIGIN = "origin"
PROJECT_ORIGIN = "PROJECT"


# Site information

INFO = """

Store and analyze metagenomic barcode data 
    """
