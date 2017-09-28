from render import metadata_loader, render_template, dict2json
import sys, os

constant_params = {
    'project' : 'metabarcode',
    'command_list' : ["setup_dir", "check_input", "create_multiqc"],
    'check_input' : 'check_input.py',
    'resdir' : './res',
    'workdir' : './work',
    'datadir' : './work/data',
    'scripts' : './scripts'

}


def setup_analysis(project,analysis_name):
    '''
    creates an analysis directory and returns its path.
    '''

    analysis_dir = project + "_" + analysis_name.strip()
    analysis_dir = os.path.join("./templates", analysis_dir)
    if not os.path.exists(analysis_dir):
        cmd = 'mkdir -p {0}'.format(analysis_dir)
        os.system(cmd)

    return analysis_dir


def get_rule(configs):
    rule = "all"
    for param in configs:
        if "action" in param.keys():
            if param["action"] == "quality report":
                rule = "qc"

        '''
        if param["action"] == "quality report":
            rule = "qc"
       '''
    return rule


def create_workflow():
    #metadata = sys.argv[1]  # metadata as a json string or json file
    #template = sys.argv[2]

    metadata = "./templates/metabarcode_qc/metabarcode_spec.json"
    template = "./templates/metabarcode_qc/metabarcode_makefile.html"
    # loads data
    configs = metadata_loader(metadata)

    parse_configs(configs)
    1/0

    rule = get_rule(configs)

    context = {'configs': configs, 'constant_params' : constant_params, 'rule': rule}

    context = dict2json(context)

    html = render_template(context, template)

    #print(html)

    # set up output directories/files.
    project = constant_params['project']
    analysis_dir = setup_analysis(project,rule)
    fname = project + "_makefile.html"
    outfile = os.path.join(analysis_dir, fname)

    with open(outfile, 'w') as fh:
        fh.write(html)


if __name__ == "__main__":
    create_workflow()
