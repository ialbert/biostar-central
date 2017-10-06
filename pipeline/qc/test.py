from pipeline import render, read_template, read_spec

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

SPEC_FILE = os.path.join(CURR_DIR, 'qc_spec.hjson' )

if __name__ =="__main__":
    spec = read_spec(SPEC_FILE)
    tmpl = read_template(spec.get("template_name", "missing"))
    html = render.render_data(spec, tmpl)
    print(html)

