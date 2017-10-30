from biostar.tools.render import render_file, read_data

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

SPEC_FILE = os.path.join(CURR_DIR, 'classify.hjson')
TEMPLATE_FILE = os.path.join(CURR_DIR, 'classify.sh')


if __name__ =="__main__":
    data = read_data(SPEC_FILE)
    html = render_file(data=data, fname=TEMPLATE_FILE)
    print(html)

