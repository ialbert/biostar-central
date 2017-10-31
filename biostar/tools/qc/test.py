
from biostar.tools.render import render_file, read_data

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))

TEMPLATE_FILE = os.path.join(CURR_DIR,"qc.sh")
SPEC_FILE = os.path.join(CURR_DIR, 'qc.hjson')


if __name__ == "__main__":
    data = read_data(SPEC_FILE)
    html = render_file(data=data, fname=TEMPLATE_FILE)
    print(html)

