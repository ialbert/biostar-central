import subprocess,sys,os


class Bunch:
    pass


def run(job):
    '''
    takes a job object, runs the job and return job status
    '''

    spec = job.spec
    template = job.template
    outdir = job.outdir
    id = job.id

    # create Makefile for analysis
    completed = subprocess.run(['python','make.py', spec, template], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    print(completed.check_returncode())

    if not os.path.exists(outdir):
            os.system('mkdir -p  {0}'.format(outdir))

    mfile = open(os.path.join(outdir, "Makefile"), "w")
    mfile.write(completed.stdout.decode('utf-8'))
    mfile.close()

    # get data


    # run make
    completed =subprocess.run(['make', 'check'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=outdir)
    print(completed.stderr.decode('utf-8'))
    print(completed.check_returncode())

    print("{0] completed successfully!".format(id))

    return


if __name__ == "__main__":

    # job creation not needed later
    spec = "templates/metabarcode_qc/metabarcode_spec.hjson"
    template = "templates/metabarcode_qc/metabarcode_makefile.html"
    out = "./workout"

    job= Bunch()
    job.id = id
    job.spec = spec
    job.template = template
    job.outdir = out

    completed =subprocess.run(['make', 'check'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=out)
    print(completed.stderr.decode('utf-8'))
    print(completed.check_returncode())
    1/0

    '''
    c=subprocess.run(['touch','hello'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=out)
    #data=c.stdout.decode('utf-8')
    #cmd="cat {0} >testfile.txt'.format(data)"
    #os.system(cmd)
    1/0
    '''

    run(job)






