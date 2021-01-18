import toml, json, os, io, base64
from django.core.management.base import BaseCommand
from django.conf import settings
from biostar.recipes.models import Analysis, Project, Data, image_path, Access
from biostar.accounts.models import User, Profile
from biostar.recipes import util, auth


class Bunch(object):
    uid = ""
    name = ""
    text = ""
    date = ""
    privacy = ""
    owner = ''
    owner_email = ''
    owner_first_name = ''
    owner_text = ''
    image = ""
    recipes = []
    json = ""
    code = ""
    is_project = False

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def parse_json(json_dict):

    uid = json_dict.get('uid')
    name = json_dict.get('name')
    text = json_dict.get('text')
    date = json_dict.get('date')
    privacy = json_dict.get('privacy')
    image = json_dict.get('image')
    owner = json_dict.get('owner_email')
    recipes = json_dict.get('recipes')
    json_text = json_dict.get('json')
    code = json_dict.get('code')
    is_project = recipes is not None

    data = Bunch(uid=uid, name=name, text=text, date=date, privacy=privacy,
                 image=image, recipes=recipes, is_project=is_project,
                 code=code, json=json_text, owner_email=owner)
    return data


def get_or_create(email, data):

    user = User.objects.get_or_create(email__iexact=email)[0]

    # Update name and text
    name = data.owner_first_name or user.profile.name
    text = data.owner_text or user.profile.text
    Profile.objects.filter(user=user).update(text=text, name=name)

    return user


def update_image(obj, img_str):

    if img_str:

        stream = io.BytesIO(base64.b64decode(img_str))

        # Convert string image to stream
        name = image_path(obj, filename="image.jpg")
        obj.image.save(name, stream, save=True)

    return


def upload_recipe(obj, project, user=None):

    recipe = Analysis.objects.filter(uid=obj.uid).first()
    json_text = toml.dumps(obj.json)

    if not recipe:
        recipe = auth.create_analysis(project=project, user=user, uid=obj.uid,
                                      name=obj.name, text=obj.text,
                                      json_text=json_text, template=obj.code)

    recipe.uid = obj.uid
    recipe.name = obj.name
    recipe.text = obj.text
    recipe.date = obj.date

    recipe.json_text = json_text or recipe.json_text
    recipe.template = obj.code
    update_image(recipe, obj.image)

    recipe.save()

    return


def upload_project(obj, user=None):

    # Check if this is a recipe or project
    project = Project.objects.filter(uid=obj.uid).first()

    if not project:
        project = auth.create_project(user=user, uid=obj.uid,
                                      name=obj.name,
                                      text=obj.text)
    project.uid = obj.uid
    project.name = obj.name
    project.text = obj.text
    project.date = obj.date
    project.privacy = obj.privacy

    update_image(project, obj.image)

    project.save()

    for recipe, vals in obj.recipes.items():
        data = parse_json(vals)
        email = data.owner_email or settings.DEFAULT_FROM_EMAIL
        owner = get_or_create(email, data=data)

        upload_recipe(data, project=project, user=owner)

    return


def upload(fname):

    stream = open(fname, 'r')

    # Get a list of projects or a single one to update.
    json_obj = json.load(stream)

    admin_email = settings.DEFAULT_FROM_EMAIL

    for key, item in json_obj.items():

        # Filter for uid if exists.
        data = parse_json(item)
        project = Project.objects.filter(uid=data.uid).first()

        email = data.owner_email or admin_email
        owner = get_or_create(email, data=data)

        if data.is_project:
            upload_project(data, user=owner)
        else:
            upload_recipe(data, project=project, user=owner)


class Command(BaseCommand):
    help = 'Interact with API end points.'

    def add_arguments(self, parser):

        # Give one file or a list of file.
        parser.add_argument('--fname', default="",
                            help="Json file to parse and upload")

        parser.add_argument('--dir', default="",
                            help="Directory of json files to upload.")

    def handle(self, *args, **options):
        fname = options['fname']
        fname = os.path.abspath(fname) if fname else ''
        directory = options['dir']
        directory = os.path.abspath(directory) if directory else ''

        if fname:
            upload(fname=fname)

        fnames = os.listdir(directory) if directory else []

        for fname in fnames:
            fname = os.path.join(directory, fname)
            upload(fname=fname)
