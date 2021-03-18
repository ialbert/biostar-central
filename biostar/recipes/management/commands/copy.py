
import logging
import psycopg2

from django.core.management.base import BaseCommand
from django.db.models import Q
from biostar.accounts.models import Profile, User
from biostar.recipes.models import Project, Data, Access, Analysis, Job
import sys


logger = logging.getLogger('engine')


def error(msg):
    logger.error(msg)
    sys.exit()


def column_names(row, colnames):
    """
    Pair items in row with corresponding column name.
    Return dictionary keyed by the column names.
    """
    mapped = {col_name: val for val, col_name in zip(row, colnames)}

    return mapped


def get_user(cursor, user_id):
    """
    Use information from source database to more_like_this users from target database.
    """
    # Get the user information from source database
    cursor.execute(f"SELECT * FROM auth_user WHERE id={user_id}")
    row = cursor.fetchone()
    colnames = [col[0] for col in cursor.description]
    obj_dict = column_names(row=row, colnames=colnames)
    user = User.objects.filter(email=obj_dict['email']).first()

    return user


def get_contributors(cursor, obj_dict):
    """
     Return owner and last edit user of object in target database
    """

    owner = get_user(user_id=obj_dict['owner_id'], cursor=cursor)
    
    if obj_dict.get('lastedit_user_id'):
        lastedit_user = get_user(user_id=obj_dict['lastedit_user_id'], cursor=cursor)
    else:
        lastedit_user = owner

    # Use the information from source database to
    # fetch users from target database.

    return owner, lastedit_user


def get_project(cursor, project_id):
    """
    Use information from source database to more_like_this projects from target database.
    """
    cursor.execute(f"SELECT * FROM engine_project WHERE id={project_id}")
    colnames = [col[0] for col in cursor.description]
    row = cursor.fetchone()
    project_dict = column_names(row=row, colnames=colnames)

    project = Project.objects.filter(uid=project_dict['uid']).first()

    return project


def get_recipe(cursor, recipe_id):
    """
    Use information from source database to more_like_this recipes from target database.
    """
    cursor.execute(f"SELECT * FROM engine_analysis WHERE id={recipe_id}")
    colnames = [col[0] for col in cursor.description]
    row = cursor.fetchone()
    recipe_dict = column_names(row=row, colnames=colnames)

    recipe = Analysis.objects.filter(uid=recipe_dict['uid']).first()

    return recipe


def copy_users(cursor):

    # Copy users over first
    users_table = 'auth_user'
    profile_table = 'accounts_profile'
    cursor.execute(f'SELECT * FROM {users_table}')
    user_rows = cursor.fetchall()
    user_colnames = [col[0] for col in cursor.description]

    logger.info(f"Copying users={len(user_rows)}.")
    for row in user_rows:

        obj_dict = column_names(row=row, colnames=user_colnames)
        user_id = obj_dict['id']
        cursor.execute(f"SELECT * FROM {profile_table} WHERE user_id={user_id}")
        # Get the first profile associated with.
        profile = cursor.fetchone()
        # Get profile obj_dict as dictionary
        profile_colnames = [col[0] for col in cursor.description]
        user_profile = column_names(row=profile, colnames=profile_colnames)

        # Update user dict with profile information.
        obj_dict.update(user_profile)
        user = User.objects.filter(Q(email=obj_dict['email']) | Q(username=obj_dict['username'])).first()
        if not user:
            # Create the user
            user = User.objects.create(username=obj_dict['username'], password=obj_dict['password'],
                                       is_superuser=obj_dict['is_superuser'], email=obj_dict['email'],
                                       is_active=obj_dict['is_active'],
                                       is_staff=obj_dict['is_staff'])
        # Update profile information
        Profile.objects.filter(user=user).update(uid=obj_dict['uid'], name=obj_dict["name"], role=obj_dict["role"],
                                                 last_login=obj_dict["last_login"], date_joined=obj_dict["date_joined"],
                                                 location=obj_dict["location"], website=obj_dict["website"],
                                                 scholar=obj_dict["scholar"], score=obj_dict["score"], text=obj_dict["text"],
                                                 twitter=obj_dict["twitter"], digest_prefs=obj_dict["digest_prefs"],
                                                 new_messages=obj_dict["new_messages"], my_tags=obj_dict["my_tags"], )
    logger.info(f"Finished copying users={len(user_rows)}.")


def copy_projects(cursor):

    projects_table = 'engine_project'
    cursor.execute(f'SELECT * FROM {projects_table}')
    project_colnames = [col[0] for col in cursor.description]
    project_rows = cursor.fetchall()

    logger.info(f"Copying projects={len(project_rows)}.")
    for row in project_rows:
        obj_dict = column_names(row=row, colnames=project_colnames)
        
        # Fetch the owner and last edit user information
        owner, lastedit_user = get_contributors(cursor=cursor, obj_dict=obj_dict)
        Project.objects.create(rank=obj_dict['rank'], lastedit_user=lastedit_user,
                               owner=owner, date=obj_dict['date'],
                               privacy=obj_dict['privacy'], uid=obj_dict['uid'], text=obj_dict['text'],
                               deleted=obj_dict['deleted'], name=obj_dict['name'],
                               lastedit_date=obj_dict['lastedit_date'], image=obj_dict['image'])

    logger.info(f"Finished copying projects={len(project_rows)}.")
    return


def copy_access(cursor):
    access_table = 'engine_access'
    cursor.execute(f'SELECT * FROM {access_table}')
    access_colnames = [col[0] for col in cursor.description]
    access_rows = cursor.fetchall()

    logger.info(f"Copying access={len(access_rows)}.")
    for row in access_rows:
        obj_dict = column_names(row=row, colnames=access_colnames)

        # Get the user and project
        user_id = obj_dict['user_id']
        project_id = obj_dict['project_id']
        user = get_user(cursor=cursor, user_id=user_id)
        project = get_project(cursor=cursor, project_id=project_id)

        access_int = obj_dict['access']
        date = obj_dict['date']

        access = Access.objects.filter(user=user, project=project).first()
        if access is None:
            Access.objects.create(user=user, project=project, access=access_int, date=date)
        else:
            Access.objects.filter(pk=access.pk).update(access=access_int, date=date)

    logger.info(f'Finished copying access={len(access_rows)}.')
    return


def copy_data(cursor):

    data_table = 'engine_data'
    cursor.execute(f'SELECT * FROM {data_table}')
    data_colnames = [col[0] for col in cursor.description]
    data_rows = cursor.fetchall()

    logger.info(f"Copying data={len(data_rows)}.")
    for row in data_rows:
        obj_dict = column_names(row=row, colnames=data_colnames)
        owner, lastedit_user = get_contributors(cursor=cursor, obj_dict=obj_dict)

        project_id = obj_dict['project_id']
        project = get_project(cursor=cursor, project_id=project_id)
        # Create the data
        Data.objects.create(method=obj_dict['method'], name=obj_dict['name'], state=obj_dict['state'],
                            image=obj_dict['image'], deleted=obj_dict['deleted'],
                            rank=obj_dict['rank'], lastedit_user=lastedit_user, owner=owner,
                            text=obj_dict['text'], date=obj_dict['date'],
                            lastedit_date=obj_dict['lastedit_date'], type=obj_dict['type'],
                            project=project, size=obj_dict['size'], file=obj_dict['file'],
                            uid=obj_dict['uid'])

    logger.info(f'Finished copying data={len(data_rows)}')

    return


def copy_analysis(cursor):

    recipe_table = 'engine_analysis'
    cursor.execute(f'SELECT * FROM {recipe_table}')
    recipe_colnames = [col[0] for col in cursor.description]
    recipe_rows = cursor.fetchall()

    logger.info(f"Copying recipes={len(recipe_rows)}.")
    for row in recipe_rows:
        obj_dict = column_names(row=row, colnames=recipe_colnames)
        owner, lastedit_user = get_contributors(cursor=cursor, obj_dict=obj_dict)
        diff_author = get_user(cursor=cursor, user_id=obj_dict['diff_author_id'])

        project_id = obj_dict['project_id']
        project = get_project(cursor=cursor, project_id=project_id)

        Analysis.objects.create(security=obj_dict['security'], deleted=obj_dict['deleted'],
                                uid=obj_dict['uid'], name=obj_dict['name'], text=obj_dict['text'],
                                owner=owner, rank=obj_dict['rank'], lastedit_user=lastedit_user,
                                lastedit_date=obj_dict['lastedit_date'],
                                diff_author=diff_author, diff_date=obj_dict['diff_date'],
                                project=project, json_text=obj_dict['json_text'],
                                template=obj_dict['template'], last_valid=obj_dict['last_valid'],
                                date=obj_dict['date'], image=obj_dict['image'])

    logger.info(f'Finished copying recipes={len(recipe_rows)}')
    return


def copy_job(cursor):
    job_table = 'engine_job'
    cursor.execute(f'SELECT * FROM {job_table}')
    job_colnames = [col[0] for col in cursor.description]
    job_rows = cursor.fetchall()

    logger.info(f"Copying jobs={len(job_rows)}.")
    for row in job_rows:
        obj_dict = column_names(row=row, colnames=job_colnames)
        owner, lastedit_user = get_contributors(cursor=cursor, obj_dict=obj_dict)
        project_id = obj_dict['project_id']
        project = get_project(cursor=cursor, project_id=project_id)

        # Get the recipe
        recipe_id = obj_dict['analysis_id']
        recipe = get_recipe(cursor=cursor, recipe_id=recipe_id)

        Job.objects.create(deleted=obj_dict['deleted'], name=obj_dict['name'], state=obj_dict['state'],
                           image=obj_dict['image'], lastedit_user=lastedit_user,
                           lastedit_date=obj_dict['lastedit_date'], owner=owner,
                           text=obj_dict['text'], date=obj_dict['date'],
                           start_date=obj_dict['start_date'], end_date=obj_dict['end_date'],
                           analysis=recipe, project=project, json_text=obj_dict['json_text'],
                           uid=obj_dict['uid'], template=obj_dict['template'],
                           security=obj_dict['security'], script=obj_dict['script'],
                           stderr_log=obj_dict['stderr_log'], stdout_log=obj_dict['stdout_log'],
                           valid=obj_dict['valid'], path=obj_dict['path'])
    logger.info(f'Finished copying jobs={len(job_rows)}')
    return


class Command(BaseCommand):
    help = "Move 'engine' table from source postgres database to current 'recipes' table."

    def add_arguments(self, parser):
        parser.add_argument('--db',  type=str, required=True, help="Source postgres database name.")
        parser.add_argument('--username', type=str, help="Username to access source postgres database")
        parser.add_argument('--password', type=str, help="Password to access source postgres database.")
        parser.add_argument('--host', type=str, help="Postgres database host.")

    def handle(self, *args, **options):
        database = options['db']
        username = options['username']
        password = options['password']
        host = options['host']

        # Connect to the database.
        try:
            conn = psycopg2.connect(database=database, user=username, password=password, host=host)
        except Exception as exc:
            logger.error(f'Error connecting to postgres database: {exc}')
            return

        cursor = conn.cursor()

        copy_users(cursor=cursor)

        copy_projects(cursor=cursor)

        copy_access(cursor=cursor)

        copy_data(cursor=cursor)

        copy_analysis(cursor=cursor)

        copy_job(cursor=cursor)
