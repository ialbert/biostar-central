from pipeline import render, read_template, read_spec

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

TEMPLATE_DIR = os.path.join(CURR_DIR,"../templates")
SPEC_FILE = os.path.join(CURR_DIR, 'qc.hjson')

if __name__ =="__main__":
    spec = read_spec(SPEC_FILE)
    #tmpl = read_template(spec.get("template_name", "missing"))
    tmpl = read_template(os.path.join(TEMPLATE_DIR, spec['template']['value']))
    html = render.render_data(spec, tmpl)
    print(html)

