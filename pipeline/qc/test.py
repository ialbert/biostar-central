from pipeline import render

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

TEMPLATE_FILE= os.path.join(CURR_DIR, 'templates', 'qc_makefile.html')
SPEC_FILE = os.path.join(CURR_DIR, 'qc_spec.hjson' )

if __name__ =="__main__":
    spec = open(SPEC_FILE).read()
    tmpl = open(TEMPLATE_FILE).read()
    html = render.render_string(spec, tmpl)
    print(html)

