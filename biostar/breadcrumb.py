from django.urls import reverse


# Constants specific to breadcrumbs builder
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


def breadcrumb_builder(icons=[], project=None, analysis=None, data=None, job=None):
    """
    This function builds the breadcrumbs on each page.
    """
    if not icons:
        return []

    path = []
    last = icons[-1]
    for icon in icons:
        is_active = icon is last

        if icon == HOME_ICON:
            step = (reverse("index"), HOME_ICON, "Home", is_active)
        elif icon == PROJECT_LIST_ICON:
            step = (reverse("project_list"), PROJECT_LIST_ICON, "Project List", is_active)
        elif icon == PROJECT_ICON:
            step = (project.url(), PROJECT_ICON, "Project View", is_active)
        elif icon == DATA_LIST_ICON:
            step = (reverse("data_list", kwargs={'uid': project.uid}), DATA_LIST_ICON, "Data Files", is_active)
        elif icon == DATA_ICON:
            step = (reverse("data_view", kwargs={'id': data.id}), DATA_ICON, f"File View", is_active)
        elif icon == DATA_UPLOAD:
            step = (reverse("data_view", kwargs={'id': project.id}), DATA_UPLOAD, f"File Upload", is_active)
        elif icon == ANALYSIS_LIST_ICON:
            step = (reverse("recipe_list", kwargs={'uid': project.uid}), ANALYSIS_LIST_ICON, "Recipe List", is_active)
        elif icon == ANALYSIS_VIEW_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_VIEW_ICON, "Recipe View", is_active)
        elif icon == ANALYSIS_RUN_ICON:
            step = (reverse("analysis_run", kwargs={'id': analysis.id}), ANALYSIS_RUN_ICON, "Analysis Run", is_active)
        elif icon == ANALYSIS_RECIPE_ICON:
            step = (reverse("recipe_view", kwargs={'id': analysis.id}), ANALYSIS_RECIPE_ICON, "Recipe Code", is_active)
        elif icon == RESULT_LIST_ICON:
            step = (reverse("job_list", kwargs={'uid': project.uid, }), RESULT_LIST_ICON, "Result List", is_active)
        elif icon == RESULT_VIEW_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_VIEW_ICON, "Result View", is_active)
        elif icon == USER_ICON:
            step = (reverse("profile"), USER_ICON, f"Profile", is_active)
        elif icon == LOGIN_ICON:
            step = (reverse("login"), LOGIN_ICON, "Login", is_active)
        elif icon == LOGOUT_ICON:
            step = (reverse("login"), LOGOUT_ICON, "Logout", is_active)
        elif icon == INFO_ICON:
            step = (reverse("info"), INFO_ICON, "Information", is_active)
        elif icon == SIGNUP_ICON:
            step = (reverse("signup"), SIGNUP_ICON, "Sign up", is_active)
        elif icon == RESULT_INDEX_ICON:
            step = (reverse("job_view", kwargs={'id': job.id}), RESULT_INDEX_ICON, "Index View", is_active)
        elif icon == ADD_USER:
            step = (reverse("project_view", kwargs={'uid': project.uid}), ADD_USER, "Manage Access", is_active)
        elif icon == PROJECT_TYPES:
            step = (reverse("project_types", kwargs={'uid': project.uid}), PROJECT_TYPES, "Manage Data Types", is_active)
        else:
            continue

        path.append(step)

    return path

