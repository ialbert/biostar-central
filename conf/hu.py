# -*- coding: utf-8 -*-
#
# import from the main settings then override some of the values
#
from main.settings import *

USE_I18N = True

USE_L10N = True

LANGUAGE_CODE = 'hu-HU'

ANON_PILL_BAR = [
    ("all", "/", "Minden", "" ),
    ("news", "/show/news/", "Hírek", "News" ),
    ("questions", "/show/questions/", "Kérdések", "Question" ),
    ("unanswered", "/show/unanswered/", "Válaszolatlan", "Unanswered" ),
    ("tutorials", "/show/tutorials/", "Lecke", "Tutorial" ),
    ("tools", "/show/tools/", "Progi", "Tool" ),
    ("videos", "/show/videos/", "Videó", "Video" ),
    ("jobs", "/show/jobs/", "Állas", "Job" ),
    ]

USER_PILL_BAR = list(ANON_PILL_BAR)

USER_PILL_BAR.insert(1,
    ("mytags", "/show/mytags/", "Saját Cimke", "" ),
)