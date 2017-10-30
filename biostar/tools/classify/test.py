from biostar.tools.render import render_data, read_data, read_template

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

SPEC_FILE = os.path.join(CURR_DIR, 'classify.hjson')
TEMPLATE_FILE = os.path.join(CURR_DIR, 'classify.sh')

if __name__ == "__main__":
    data = read_data(fname=SPEC_FILE)
    template = read_template(fname=TEMPLATE_FILE)
    html = render_data(context=data, template_txt=template)
    print(html)


