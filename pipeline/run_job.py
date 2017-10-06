import subprocess, os
import sys, json, hjson
from pipeline.render import render_file,render_string

CURR_DIR = os.path.dirname(os.path.realpath(__file__))


def run(job):
    ''''
    takes job object, runs the job and return job status
    '''

    spec = job.spec
    template = job.template
    outdir = job.outdir
    #jobid = job.jobid

    errorlog = []

    try:

        mtext = render_string(spec,template)

        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        with open(os.path.join(outdir, "Makefile"), 'wt') as fp:
            fp.write(mtext)

        process = subprocess.run(['make', 'all'], cwd=outdir, stderr=subprocess.PIPE, check=True)
        job.status = process.returncode

    except subprocess.CalledProcessError as err:
        print("ERROR!!")
        errorlog.append(err.stderr.decode('utf-8'))
        print("***",err.returncode)
        job.status = err.returncode

    finally:
        print(CURR_DIR)
        #errorlog.append(process.stderr.decode('utf-8'))
        job.log = "\n".join(errorlog)
        print("printing errlog")
        print(job.log)
        #job.save()
        #print("Job Done")

    return job


def get_spec(specfile='spec'):
    spec = hjson.load(open(specfile))
    return spec


def test_job(data='data', sinfo='sinfo'):

    spec_file = "templates/qc/qc_spec.hjson"
    template_file = "templates/qc/qc_makefile.html"


    spec = get_spec(specfile=spec_file)

    spec['data']['value'] = data
    spec['sampleinfo']['value'] = sinfo

    new_spec = hjson.dumps(spec)
    template_txt = open(template_file).read()

    #html=render_string(new_spec,template_txt)
    #print(html)
    #1/0

    job = Job()
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

    test_job(data="~/work/web-dev/biostar-engine/pipeline/testdata/data.tar.gz",sinfo="~/work/web-dev/biostar-engine/pipeline/testdata/sampleinfo.txt")
    #run(job)






