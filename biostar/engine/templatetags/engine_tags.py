from textwrap import dedent
import json
import logging
from django import template
from django.contrib.staticfiles.templatetags.staticfiles import static
from django.template import defaultfilters
from django.utils.safestring import mark_safe
from django.urls import reverse
from django.db.models import Q
from django.test.client import RequestFactory

from biostar import settings
from biostar.engine.models import Job, make_html, Project, Data, Analysis, Access
from biostar.engine import auth


logger = logging.getLogger("engine")
register = template.Library()

JOB_COLORS = {
    Job.SPOOLED: "violet",
    Job.ERROR: "red", Job.QUEUED: "teal",
    Job.RUNNING: "orange", Job.COMPLETED: "green"
}

DATA_COLORS = {
    Data.PENDING: "teal", Data.READY:"green", Data.ERROR :"red"
}



@register.simple_tag
def moderate(request):

    url = reverse("recipe_mod")
    user = request.user
    recipes = Analysis.objects.filter(security=Analysis.UNDER_REVIEW).count()

    correct = f"<b>{recipes} recipes</b> need" if recipes > 1 else f"<b>{recipes} recipe</b> needs"
    template = ''
    if user.is_authenticated and user.profile.is_moderator and (recipes > 0):

        template = f"""
            <div class="ui mini message">
            <i class="ui check circle icon"></i>
            {correct} to be <a href={url}>reviewed and authorized</a>.
            </div>
            """

    return mark_safe(template)


@register.simple_tag
def sticky_label(obj):
    label = mark_safe('<span class ="ui label">Sticky</span>')
    return label if obj.sticky else ''

@register.simple_tag
def build_path(path, name):
    return f'{path}/{name}'

@register.inclusion_tag('widgets/job_file_list.html', takes_context=True)
def job_file_list(context, path, files, job, form=None):
    back = "/".join(path.split("/")[:-1])
    return dict(path=path, files=files, job=job, form=form, back=back)

@register.inclusion_tag('widgets/file_list.html', takes_context=True)
def file_list(context, path, files, obj, form=None):

    back = "/".join(path.split("/")[:-1])

    if isinstance(obj, Data):
        view_url, serve_url = 'data_view', 'data_serve'
    else:
        view_url, serve_url = 'job_view', 'job_serve'

    return dict(path=path, files=files, obj=obj, form=form, back=back, view_url=view_url, serve_url=serve_url)


@register.inclusion_tag('widgets/action_bar.html', takes_context=True)
def action_bar(context, instance, edit_url, extra=""):

    edit_url = reverse(edit_url, kwargs=dict(uid=instance.uid))

    if isinstance(instance, Job):
        obj_type = "job"
    elif isinstance(instance, Data):
        obj_type = "data"
    elif isinstance(instance, Analysis):
        obj_type = "recipe"
    else:
        obj_type = None

    action_url = reverse("toggle_state", kwargs=dict(uid=instance.uid, obj_type=obj_type))
    return dict(instance=instance, edit_url=edit_url, user=context["user"],
                action_url=action_url, extra=extra)


@register.inclusion_tag('widgets/list_view.html')
def list_view(projects=None, data_list=None, recipe_list=None, job_list=None):

    return dict(projects=projects, data_list=data_list, recipe_list=recipe_list,
                job_list=job_list)


@register.inclusion_tag('widgets/recipe_moderate.html')
def recipes_moderate(cutoff=0):

    recipes = Analysis.objects.filter(security=Analysis.UNDER_REVIEW,
                                      deleted=False).order_by("-date")

    cutoff = cutoff or len(recipes)
    return dict(recipes=recipes[:cutoff])


@register.filter
def is_checkbox(field):
    "Check if current field is a checkbox"

    return True if field.field.widget.input_type == "checkbox" else False


@register.filter
def has_files(request):
    "Checks if object in clipboard is a list of files belonging to a job"
    files = request.session.get("files_clipboard", None)

    # Last item loaded into files clipboard is the instance.uid the files belong to.
    root_path = auth.validate_files_clipboard(request=request)
    if not root_path:
        return False

    # Some files in clipboard might be outside job path.
    files = [f for f in files if f.startswith(root_path)]

    # Files might still be empty at this point
    return True if files else False


@register.simple_tag
def get_projects(user):

    "Used to return projects list in the profile."

    projects = auth.get_project_list(user=user, include_public=False)


    access = Access.objects.filter(project__in=projects, user=user,
                                   access__gt=Access.NO_ACCESS)

    return [(access, access.project) for access in access]


def update_dict(iter):

    results=[dict(
                title=d.name,
                description=d.summary,
                url=d.url()) for d in iter
            ]

    return results


@register.inclusion_tag('widgets/search.html')
def search(request):

    projects = auth.get_project_list(user=request.user)

    data = Data.objects.filter(deleted=False,project__in=projects)
    recipe = Analysis.objects.filter(deleted=False, project__in=projects)
    jobs = Job.objects.filter(deleted=False, project__in=projects)
    projects = projects.filter(deleted=False)

    content = [dict(name="Projects",results =update_dict(iter=projects)),
               dict(name="Data", results=update_dict(iter=data)),
               dict(name="Recipes", results=update_dict(iter=recipe)),
               dict(name="Jobs", results=update_dict(iter=jobs)),
               ]

    return dict(content=json.dumps(content))


@register.inclusion_tag('widgets/paste.html')
def paste(action_view, obj, form, contains="Nothing"):
    "Default provides template for pasting a recipe"

    action_url = reverse(action_view, kwargs=dict(uid=obj.uid))


    return dict(action_url=action_url, form=form, contains=contains)


@register.filter
def has_recipe(request):
    """
    Checks if object in clipboard is a recipe.
    """
    uid = request.session.get("recipe_clipboard")
    recipe = Analysis.objects.filter(uid=uid).first()
    return bool(uid and recipe)



@register.simple_tag
def privacy_label(project):
    label = mark_safe(f'<span class ="ui label">{project.get_privacy_display()}</span>' )
    return label


@register.inclusion_tag('widgets/authorization_required.html')
def security_label(analysis):
    context = dict(analysis=analysis)
    return context

@register.simple_tag
def job_color(job):
    """
    Returns a color based on job status.
    """
    return JOB_COLORS.get(job.state, "")

@register.simple_tag
def activate(value1, value2):
    """
    Returns a color based on job status.
    """
    return "active" if value1  == value2 else ''


@register.simple_tag
def data_color(data):
    "Return a color based on data status."

    return DATA_COLORS.get(data.state,  "")

@register.simple_tag
def access_color(user, project):

    if user.is_authenticated:
        access = Access.objects.filter(user=user, project=project).first()
    else:
        access = None

    if access and access.access in (Access.WRITE_ACCESS, Access.OWNER_ACCESS):
        return "green"
    else:
        return ""


@register.simple_tag
def type_label(data):

    if data.type:
        label = lambda x:f"<span class='ui label' > {x} </span>"
        types = [label(t) for t in data.type.split(',')]
        return mark_safe(''.join(types))
    return ""

@register.simple_tag
def state_label(data, error_only=False):

    label = f'<span class="ui { DATA_COLORS.get(data.state, "") } label"> {data.get_state_display()} </span>'

    # Error produce error only.
    if error_only and data.state != Data.ERROR:
        label = ""

    return mark_safe(label)


@register.simple_tag
def img(obj):
    """
    Returns the image associated with the object or a placeholder
    """
    if obj.image:
        return obj.image.url
    else:
        return static("images/placeholder.png")


@register.inclusion_tag('widgets/show_messages.html')
def show_messages(messages):
    """
    Renders the messages
    """
    return dict(messages=messages)


@register.inclusion_tag('widgets/project_name_bar.html', takes_context=True)
def project_name_bar(context, project):
    """
    Returns a label for project.
    """
    return dict(project=project)


@register.inclusion_tag('widgets/recipe_form.html')
def recipe_form(form):
    """
    Renders a recipe form.
    """
    return dict(form=form)

@register.inclusion_tag('widgets/access_form.html')
def access_form(project, user, form):
    """
    Generates an access form.
    """

    return dict(project=project, user=user, form=form)


@register.inclusion_tag('widgets/project_action_bar.html')
def project_action_bar(user, project):
    return dict(use=user, project=project)

@register.inclusion_tag('widgets/job_elapsed.html')
def job_minutes(job):

    return dict(job=job)


@register.simple_tag
def size_label(data):
    """
    Returns a label for data sizes.
    """
    size = f"{defaultfilters.filesizeformat(data.size)}"
    return mark_safe(f"<span class='ui mini label'>{size}</span>")


@register.inclusion_tag('widgets/form_errors.html')
def form_errors(form):
    """
    Turns form errors into a data structure
    """
    try:
        errorlist = [ ('', message) for message in form.non_field_errors() ]

        for field in form:
            for error in field.errors:
                errorlist.append((f'{field.name}:', error))
    except Exception:
        errorlist = []

    context = dict(errorlist=errorlist)

    return context


@register.filter
def markdown(text):
    """
    Generates HTML from a markdown value.
    """
    if not text:
        return ''
    text = dedent(text)
    html = make_html(text)
    return mark_safe(html)


@register.inclusion_tag('widgets/menubar.html', takes_context=True)
def menubar(context, request=None):
    user = context.request.user

    return dict(user=user, request=request)
