from fabric.api import *

main_path = '/export/sites/main_engine'
test_path = '/export/sites/test_engine'

env_prefix = ('source ~/miniconda3/envs/engine/bin/activate engine')

def remote_pull(path):
    with cd(path), prefix(env_prefix):
        run('git pull')

def test_pull():
    remote_pull(test_path)

def main_pull():
    remote_pull(main_path)

def restart_nginx():
    sudo("service nginx restart")

def restart_uwsgi():
    sudo("supervisorctl restart all")