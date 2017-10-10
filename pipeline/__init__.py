import os
import hjson

__THIS = os.path.dirname(__file__)

def read_template(path):
    fullpath = os.path.join(__THIS, path)
    content = open(fullpath).read()
    return content

def read_spec(path):
    data = hjson.load(open(path))
    return data
