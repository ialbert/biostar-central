from pipeline import render, read_template, read_spec

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

SPEC_FILE = os.path.join(CURR_DIR, 'classify.hjson')
TEMPLATE_FILE = os.path.join(CURR_DIR, 'classify_makefile.txt')


if __name__ =="__main__":
    spec = read_spec(SPEC_FILE)
    #tmpl = read_template(spec.get("template_name", "missing"))
    tmpl = read_template(TEMPLATE_FILE)
    html = render.render_data(spec, tmpl)
    print(html)

