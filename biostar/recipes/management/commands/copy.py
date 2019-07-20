
import logging
import psycopg2

from django.core.management.base import BaseCommand
from django.db.models import Q
from biostar.accounts.models import Profile, User
from biostar.recipes.models import Project
import sys


logger = logging.getLogger('engine')


def error(msg):
    logger.error(msg)
    sys.exit()


def column_names(row, colnames):
    """
    Pair items in row with corresponding column name.
    Return dictionary keyed by column names.
    """
    mapped = {col_name: row for row, col_name in zip(row, colnames)}

    return mapped


def get_contributors(cursor, obj_dict):
    """
     Return owner and last edit user of object in target database
    """
    def get_user(user_id):
        # Get the user information from source database
        cursor.execute(f"SELECT * FROM auth_user WHERE id={user_id}")
        row = cursor.fetchone()
        colnames = [col[0] for col in cursor.description]
        obj_dict = column_names(row=row, colnames=colnames)
        return obj_dict
    
    owner_dict = get_user(user_id=obj_dict['owner_id'])
    
    if obj_dict.get('lastedit_user_id'):
        lastedit_user_dict = get_user(user_id=obj_dict['lastedit_user_id'])
    else:
        lastedit_user_dict = owner_dict

    # Use the information from source database to
    # fetch users from target database.
    lastedit_user = User.objects.filter(email=lastedit_user_dict['email']).first()
    owner = User.objects.filter(email=owner_dict['email']).first()

    return owner, lastedit_user


def copy_users(cursor):

    # Copy users over first
    users_table = 'auth_user'
    profile_table = 'accounts_profile'
    cursor.execute(f'SELECT * FROM {users_table}')
    user_rows = cursor.fetchall()
    user_colnames = [col[0] for col in cursor.description]

    logger.info(f"Copying {len(user_rows)} users.")
    for row in user_rows:
        obj_dict = column_names(row=row, colnames=user_colnames)
        user_id = obj_dict['id']
        cursor.execute(f"SELECT * FROM {profile_table} WHERE user_id={user_id}")
        # GEt the first profile associated with.
        profile = cursor.fetchone()
        # Get profile obj_dict as dictionary
        profile_colnames = [col[0] for col in cursor.description]
        user_profile = column_names(row=profile, colnames=profile_colnames)

        # Update user dict with profile obj_dict.
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
    logger.info(f"Finished copying {len(user_rows)} users.")


def copy_projects(cursor):

    projects_table = 'engine_project'
    cursor.execute(f'SELECT * FROM {projects_table}')
    project_colnames = [col[0] for col in cursor.description]
    project_rows = cursor.fetchall()

    logger.info("Copying all projects.")
    for row in project_rows:
        obj_dict = column_names(row=row, colnames=project_colnames)
        
        # Fetch the owner and last edit user information
        owner, lastedit_user = get_contributors(cursor=cursor, obj_dict=obj_dict)
        Project.objects.create(rank=obj_dict['rank'], lastedit_user=lastedit_user,
                               owner=owner, date=obj_dict['date'],
                               privacy=obj_dict['privacy'], uid=obj_dict['uid'], text=obj_dict['text'],
                               deleted=obj_dict['deleted'], name=obj_dict['name'],
                               lastedit_date=obj_dict['lastedit_date'], image=obj_dict['image'])

    logger.info(f"Finished copying {len(project_rows)} projects.")
    return


def copy_data(cursor):

    data_table = 'engine_data'
    cursor.execute(f'SELECT * FROM {data_table}')
    data_colnames = [col[0] for col in cursor.description]
    data_rows = cursor.fetchall()



    return


def copy_analysis(cursor):
    return


def copy_job():
    return


def copy_access():
    return


class Command(BaseCommand):
    help = "Move 'engine' table from source postgres database to current 'recipes' table."

    def add_arguments(self, parser):
        parser.add_argument('--db',  type=str, required=True, help="Source postgres database name.")
        parser.add_argument('--username', type=str, help="Username to access source postgres database")
        parser.add_argument('--password', type=str, help="Password to access source postgres database.")

    def handle(self, *args, **options):
        database = options['db']

        # Connect to the database.

        try:
            conn = psycopg2.connect(database=database)
        except Exception as exc:
            logger.error(f'Error connecting to postgres database: {exc}')
            return

        cursor = conn.cursor()

        copy_users(cursor=cursor)

        copy_projects(cursor=cursor)

        # Select every table belonging to 'engine' app

        # Add values to current databases 'recipes' table
