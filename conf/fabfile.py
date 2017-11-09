from fabric.api import *
import sys


if tuple(sys.version_info) > (2, 8):
    print("Error: Fabric requires Python 2.7")
    sys.exit(-1)

main_path = '/export/sites/main_engine'
test_path = '/export/sites/test_engine'
psu_path = '/export/sites/psu_engine'

test_env = ('source ~/miniconda3/envs/engine/bin/activate test_env')
main_env = ('source ~/miniconda3/envs/engine/bin/activate main_env')
psu_env = ('source ~/miniconda3/envs/psu_env/bin/activate psu_env')



def deploy_latest(path, env, name='main'):
    with cd(path), prefix(env):
        run('git pull')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        sudo("supervisorctl restart %s" % name)

def deploy_reset(path, env, name='test'):
    with cd(path), prefix(env):
        run("rm -f export/engine.db")
        run('git pull')
        #run('make install')
        run('make clean testdata')
        run('make reset')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        sudo('supervisorctl restart %s' % name)


def deploy_reset_test():
    deploy_reset(test_path, test_env, name='test')


def deploy_reset_main():
    deploy_reset(main_path, main_env, name='main')

def deploy_test():
    deploy_reset(test_path, test_env, name='test')

def deploy_main():
    deploy_reset(main_path, main_env, name='main')

def deploy_psu():
    deploy_reset(psu_path, psu_env, name='main')


def restart_nginx():
    sudo("service nginx restart")

def restart_all():
    sudo("supervisorctl restart all")
