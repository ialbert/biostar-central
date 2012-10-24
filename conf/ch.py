# -*- coding: utf-8 -*-
#
# import from the main settings then override some of the values
#A
from main.settings import *

USE_I18N = True

USE_L10N = True

LANGUAGE_CODE = 'zh-CH'

# to override a location put it before the default templates
# usually needs full paths
TEMPLATE_DIRS = [ 'conf/custom-html' ] + TEMPLATE_DIRS


ANON_PILL_BAR = [
    ("all", "/", "所有", "" ),
    ("news", "/show/news/", "新闻", "News" ),
    ("questions", "/show/questions/", "问题", "Question" ),
    ("unanswered", "/show/unanswered/", "未解答", "Unanswered" ),
    ("tutorials", "/show/tutorials/", "课程", "Tutorial" ),
    ("tools", "/show/tools/", "工具", "Tool" ),
    ("videos", "/show/videos/", "视频", "Video" ),
    ("jobs", "/show/jobs/", "招聘", "Job" ),
    ]

USER_PILL_BAR = list(ANON_PILL_BAR)

USER_PILL_BAR.insert(1,
    ("mytags", "/show/mytags/", "我的标签", "" ),
)