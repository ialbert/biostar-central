from django.urls import reverse


# Constants specific to breadcrumb
HOME_ICON = "home"
PROJECT_LIST_ICON = "database"
PROJECT_ICON = "archive"
PROJECT_TYPES = "file outline"
DATA_LIST_ICON = "file text"
DATA_ICON = "file"
ANALYSIS_LIST_ICON = "settings"
ANALYSIS_VIEW_ICON = "setting"
ANALYSIS_RUN_ICON = "spinner"
ANALYSIS_RECIPE_ICON = "list layout"
RESULT_LIST_ICON = "bar chart"
RESULT_VIEW_ICON = "area chart"
RESULT_INDEX_ICON = "list layout icon"
LOGIN_ICON = "sign in"
LOGOUT_ICON = "sign out"
INFO_ICON = "info circle icon"
SIGNUP_ICON = "add user icon"
USER_ICON = "user icon"
DATA_UPLOAD="upload icon"
ADD_USER = "users icon"


class Bunch(object):

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


# Placeholder to mimic attrs
T = Bunch(url=lambda:None, id=0, uid="0")


def icon_mapper(project=None, analysis=None, job=None, data=None):
    "Map an icon to its info ( url and name)"

    return {    HOME_ICON: dict(url=reverse("index"), name="Home"),
                PROJECT_LIST_ICON: dict(url=reverse("project_list"), name="Project List"),
                PROJECT_ICON:dict(url=project.url(), name="Project View"),
                ADD_USER: dict(url=project.url(), name="Manage Access"),
                USER_ICON: dict(url=reverse("profile"), name="Profile"),
                LOGIN_ICON: dict(url=reverse("login"), name="Login"),
                LOGOUT_ICON: dict(url=reverse("logout"), name="Logout"),
                SIGNUP_ICON: dict(url=reverse("signup"), name="Sign up"),
                DATA_LIST_ICON: dict(url=reverse("data_list", kwargs={'uid': project.uid}), name="Data Files"),
                DATA_ICON: dict(url=reverse("data_view", kwargs={'id': data.id}), name="File View"),
                DATA_UPLOAD: dict(url=reverse("data_view", kwargs={'id': data.id}), name="File Upload"),
                ANALYSIS_LIST_ICON: dict(url=reverse("recipe_list", kwargs={'uid': project.uid}), name="Recipe List"),
                ANALYSIS_VIEW_ICON: dict(url=reverse("recipe_view", kwargs={'id': analysis.id}), name="Recipe View"),
                ANALYSIS_RUN_ICON: dict(url=reverse("analysis_run", kwargs={'id': analysis.id}), name="Analysis Run"),
                ANALYSIS_RECIPE_ICON: dict(url=reverse("recipe_view", kwargs={'id': analysis.id}), name="Recipe Code"),
                RESULT_LIST_ICON: dict(url=reverse("job_list", kwargs={'uid': project.uid}), name="Results List"),
                RESULT_VIEW_ICON: dict(url=reverse("job_view", kwargs={'id': job.id}), name="Result View"),
                RESULT_INDEX_ICON: dict(url=reverse("job_view", kwargs={'id': job.id}), name="Index View"),
                PROJECT_TYPES: dict(url=reverse("project_types", kwargs={'uid': project.uid}),
                                    name="Manage Data Types")
                }

def breadcrumb_builder(icons=[], project=T, analysis=T, job=T, data=T):
    """
    This function builds the breadcrumbs on each page.
    """
    if not icons:
        return []

    path = []
    last = icons[-1]

    for icon in icons:

        is_active = icon is last

        ready = icon_mapper(project=project, analysis=analysis, job=job, data=data).get(icon)
        if not ready:
            continue

        step = (ready.get("url", "/"), icon, ready.get("name", "NOT FOUND"), is_active)
        path.append(step)

    return path
