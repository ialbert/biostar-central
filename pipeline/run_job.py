import subprocess,os
import urllib.request
#from engine.models import Job


DATA="http://iris.bx.psu.edu/projects/metabarcode-data/data.tar.gz"
SINFO="http://iris.bx.psu.edu/projects/metabarcode-data/sampleinfo.txt"


def run(job):
    '''
    takes a job object, runs the job and return job status
    '''

    spec = job.spec
    template = job.template
    outdir = job.outdir
    jobid = job.jobid

    # create Makefile for analysis
    process1 = subprocess.run(['python','make.py', spec,template], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(process1.check_returncode())
    
    if not os.path.exists(outdir):
            os.system('mkdir -p  {0}'.format(outdir))

    mfile = open(os.path.join(outdir, "Makefile"), "w")
    mfile.write(process1.stdout.decode('utf-8'))
    mfile.close()

    # get analysis data
    urllib.request.urlretrieve(SINFO, os.path.join(outdir, 'sampleinfo.txt'))
    urllib.request.urlretrieve(DATA, os.path.join(outdir, 'data.tar.gz'))
    
    # run analysis
    process2 = subprocess.run(['make', 'all'], cwd=outdir)
    print(process2.check_returncode())
    
    print("{0} ran successfully!".format(jobid))
    return


class Job:
    pass


if __name__ == "__main__":

    # these won't be needed later.
    ############# 
    spec = "templates/metabarcode_qc/metabarcode_spec.hjson"
    template = "templates/metabarcode_qc/metabarcode_makefile.html"
    out = "./workout"
    jobid = "job0"

    job = Job()
    job.jobid= jobid
    job.spec = spec
    job.template = template
    job.outdir = out
    ##########

    run(job)






