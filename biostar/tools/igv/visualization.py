from biostar.tools.render import render_template
import sys, os, glob


def make_igv(url, files_dir):

    # Get the genome file.
    fasta_idx = os.path.join(files_dir, f'*.fa.fai')
    genome = glob.glob(fasta_idx)[0].replace(".fai","")
    genome = os.path.basename(genome)

    # Collect all bam files.
    bam_patt = os.path.join(files_dir, f'*.bam')
    bam_files = glob.glob(bam_patt)

    def get_fname(file_list):
        fnames = [ os.path.splitext(os.path.basename(f))[0] for f in file_list]
        return fnames

    samples = get_fname(bam_files)

    path = os.path.join(url,os.path.basename(files_dir))

    igv_template = "igv.xml"
    data = dict(genome=genome, path=path, samples=samples)
    html = render_template(data, igv_template)
    print(html)


if __name__ == "__main__":
    url = sys.argv[1]
    workdir = sys.argv[2]
    make_igv(url, workdir)