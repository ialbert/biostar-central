#
# import from the main settings then override some of the values
#
from main.settings import *

EXTERNAL_AUTHENICATION = {
    "TEST-KEY" : ("abcd", "User <b>%(display_name)s</b> with email <b>%(email)s</b> has the tags <b>%(tags)s</b>"),
}