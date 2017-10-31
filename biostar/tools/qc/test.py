from biostar.tools.render import read_data, read_template, render_data

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

TEMPLATE_FILE = os.path.join(CURR_DIR, "qc.sh")
SPEC_FILE = os.path.join(CURR_DIR,  'qc.hjson')


if __name__ == "__main__":
    data = read_data(fname=SPEC_FILE)
    template = read_template(fname=TEMPLATE_FILE)
    html = render_data(context=data, template_txt=template)
    print(html)

