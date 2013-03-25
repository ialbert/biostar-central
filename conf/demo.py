#
# import from the main settings then override some of the values
#
from main.settings import *

patt = """Galaxy username <b>%(display_name)s</b>\n
Tool name: <b>%(tool_name)s</b>
Tool version:  <b>%(tool_version)s</b>
Tool id: <b>%(tool_id)s</b>
Tags: <b>%(tags)s</b>"""

EXTERNAL_AUTHENICATION = {
    "TEST-KEY" : ("abcd", patt),
}