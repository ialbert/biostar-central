import subprocess, os
import sys, json, hjson
#from pipeline.render import render_file,render_string
from pipeline import render, read_template, read_spec

import os

CURR_DIR = os.path.dirname(os.path.realpath(__file__))
SPEC_FILE = os.path.join(CURR_DIR, 'qc_spec.hjson' )


def run(job):
    ''''
    takes job object, runs the job and return job status
    '''

    spec = job.spec
    template = job.template
    outdir = job.outdir
    jobid = job.jobid

    errorlog = []

    try:

        mtext = render.render_data(spec,template)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        with open(os.path.join(outdir, "Makefile"), 'wt') as fp:
            fp.write(mtext)

        process = subprocess.run(['make', 'all'], cwd=outdir, stderr=subprocess.PIPE, check=True)
        job.status = process.returncode

    except subprocess.CalledProcessError as err:
        print("ERROR!!")
        errorlog.append(err.stderr.decode('utf-8'))
        if len(errorlog) > 100:
            sys.exit(1)
        job.status = err.returncode

    finally:
        print(CURR_DIR)
        job.log = "\n".join(errorlog)
        print("printing errlog")
        print(job.log)
        #job.save()
    return job


    job.jobid= "job0"
    job.spec = new_spec
    job.template = template_txt
    job.outdir = "./work"
    job.status =""


    run(job)

    print (job.status)
    print (job.log)


class Job:
    pass


if __name__ == "__main__":
    spec = read_spec(SPEC_FILE)
    tmpl = read_template(spec.get("template_name", "missing"))

    job = Job()
    job.spec  = spec
    job.template = tmpl
    job.outdir = "./"
    job.jobid= "job0" 	
    job.status = None
    
    run(job)



