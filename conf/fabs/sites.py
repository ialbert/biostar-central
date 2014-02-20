

from fabric.api import *

def usegalaxy(user="www"):
    "Sets the environment for the biostar galaxy"
    env.hosts = ['biostar.usegalaxy.org']
    env.user = user

def hostname():
    run("hostname")
