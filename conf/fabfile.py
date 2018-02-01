from fabric.api import *
import sys

if tuple(sys.version_info) > (2, 8):
    print("Error: Fabric requires Python 2.7")
    sys.exit(-1)

# The path to the engine.
engine_path = '/export/sites/biostar-engine'

# The path to the recipes.
recipe_path = '/export/sites/biostar-recipes'

# The default conda environment.
conda_env = ('source ~/miniconda3/envs/engine/bin/activate engine')

def restart_nginx():
    sudo("service nginx restart")

def restart_server():
    sudo("supervisorctl reload")

def update_recipes(path):
    """
    Updates the recipe directory
    """
    with cd(recipe_path):
        run('git pull')

    with cd(path):
        run ('make data')

def deploy(path=engine_path, env=conda_env, name='engine', recipes=False):

    # Should we reset the engine.
    if recipes:
        update_recipes(path=path)

    with cd(path), prefix(env):
        sudo('supervisorctl stop %s' % name)
        run('git pull')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        sudo("supervisorctl restart %s" % name)

def sqlite_reset(path=engine_path, env=conda_env, name='engine'):

    with cd(path), prefix(env):
        sudo('supervisorctl stop %s' % name)
        run("make reset")

    deploy(path=path, env=env, name=name)



