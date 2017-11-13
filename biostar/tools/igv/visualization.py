from biostar.tools.render import render_template
import sys, os, glob


FILE_EXTENSIONS=[".bam", ".bw", ".vcf", ".gff3", ".gtf", ".bedgraph", "bigwig"]


def get_files(files_dir, file_type):
    '''
    returns a list of files of specific type
    '''
    file_patt = os.path.join(files_dir, file_type)
    selected_files = glob.glob(file_patt)
    return selected_files


def make_igv(url, files_dir):

    # Get all file types in files_dir.
    file_list = list(set(os.path.splitext(file)[1] for file in os.listdir(files_dir)))
    ext_list= [x for x in file_list if x in FILE_EXTENSIONS]

    # Get track types
    track_types = list(map(lambda item: item.replace(".", ""), ext_list))

    # Get the genome file.
    fasta_idx = os.path.join(files_dir, f'*.fa.fai')
    genome = glob.glob(fasta_idx)[0].replace(".fai","")
    genome = os.path.basename(genome)

    # Collect all bam files.
    bam_files = get_files(files_dir,f'*.bam')

    # Collect all samples.
    samples = [ os.path.splitext(os.path.basename(f))[0] for f in bam_files]

    path = os.path.join(url,os.path.basename(files_dir))

    # Render template.
    igv_template = "igv.xml"
    data = dict(genome=genome, path=path, samples=samples, track_types=track_types)
    html = render_template(data, igv_template)
    print(html)


if __name__ == "__main__":
    url = sys.argv[1]
    workdir = sys.argv[2]
    make_igv(url, workdir)