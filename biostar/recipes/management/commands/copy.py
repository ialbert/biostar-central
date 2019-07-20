
import logging
import psycopg2

from django.core.management.base import BaseCommand

from biostar.accounts.models import Profile, User
import sys


logger = logging.getLogger('engine')


def error(msg):
    logger.error(msg)
    sys.exit()


#def colnames():
#    return


def column_names(row, colnames):
    """
    Pair items in row with corresponding column name.
    Return dictionary keyed by column names.
    """
    mapped = {col_name: row for row, col_name in zip(row, colnames)}

    return mapped


def copy_users(cursor):

    # Copy users over first
    users_table = 'auth_user'
    profile_table = 'accounts_profile'
    cursor.execute(f'SELECT * FROM {users_table}')
    user_cols = [col[0] for col in cursor.description]

    user_rows = cursor.fetchall()

    logger.info("Copying all users.")
    for row in user_rows:
        user = column_names(row=row, colnames=user_cols)
        cursor.execute(f"SELECT * FROM {profile_table} WHERE user_id={user['id']}")

        profile = cursor.fetchone()
        # Get column names for profile table

        profile_colnames = [col[0] for col in cursor.description]
        # Get profile info as dictionary
        user_profile = column_names(row=profile, colnames=profile_colnames)
        # Update user dict with profile info.
        user.update(user_profile)

        added_user = User.objects.create(username=user['username'], password=user['password'],
                                         is_superuser=user['is_superuser'], email=user['email'],
                                         is_active=user['is_active'],
                                         is_staff=user['is_staff'])

        # Update user profile
        Profile.objects.filter(user=added_user).update(uid=user["uid"], name=user["name"], role=user["role"],
                                                       last_login=user["last_login"], date_joined=user["date_joined"],
                                                       location=user["location"], website=user["website"],
                                                       scholar=user["scholar"], score=user["score"], text=user["text"],
                                                       twitter=user["twitter"], digest_prefs=user["digest_prefs"],
                                                       new_messages=user["new_messages"], my_tags=user["my_tags"], )
    logger.info("Finished copying all users.")


def copy_projects(cursor):

    projects_table = 'engine_project'
    cursor.execute(f'SELECT * FROM {projects_table}')

    project_cols = [col[0] for col in cursor.description]
    project_rows = cursor.fetchall()

    for row in project_rows:
        project = column_names(row=row, colnames=project_cols)
        print(project)

    return


class Command(BaseCommand):
    help = "Move 'engine' table from source postgres database to current 'recipes' table."

    def add_arguments(self, parser):
        parser.add_argument('--db',  type=str, required=True, help="Existing postgres database name.")
        parser.add_argument('--username', type=str, help="Username to access postgres database")
        parser.add_argument('--password', type=str, help="Password to access postgres database.")

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

        #copy_projects(cursor=cursor)

        # Select every table belonging to 'engine' app

        # Add values to current databases 'recipes' table
