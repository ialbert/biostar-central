from fabric.api import *
import sys

if tuple(sys.version_info) > (2, 8):
    print("Error: Fabric requires Python 2.7")
    sys.exit(-1)

main_path = '/export/sites/main_web'
test_path = '/export/sites/test_web'
psu_path = '/export/sites/psu_engine'

test_env = ('source ~/miniconda3/envs/test_env/bin/activate test_env')
main_env = ('source ~/miniconda3/envs/main_env/bin/activate main_env')
psu_env = ('source ~/miniconda3/envs/psu_env/bin/activate psu_env')

def restart_nginx():
    sudo("service nginx restart")

def restart_all():
    sudo("supervisorctl restart all")

def deploy(path, env, name='main'):
    with cd(path), prefix(env):
        run('git pull')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        sudo("supervisorctl restart %s" % name)

def reset(path, env, name='test'):
    with cd(path), prefix(env):
        sudo('supervisorctl stop %s' % name)
        run("rm -f export/engine.db")
        run('git pull')
        #run ('make conda')
        #run('make install')
        #run('make testdata')
        run('make reset')
        run('python manage.py migrate')
        run('python manage.py collectstatic --noinput -v 0')
        sudo('supervisorctl start %s' % name)

def deploy_main():
    deploy(main_path, main_env, name='main')

def deploy_psu():
    deploy(psu_path, psu_env, name='main')
    #reset(psu_path, psu_env, name='main')


