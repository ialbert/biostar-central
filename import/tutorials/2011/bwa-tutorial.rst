Title: Mapping reads with bwa and bowtie
Tags: bowtie bwa tutorial MSU-NGS-2011
##rest

"Mapping reads with bwa and bowtie tutorial" imported from the MSU course **Analyzing Next Generation Sequencing Data** (http://bioinformatics.msu.edu/ngs-summer-course-2011)

In this tutorial, we're going to take a set of Illumina reads from an inbred Drosophila melanogaster line, 
and map them back to the reference genome. (After these steps, we could do things like generate a list of SNPs 
at which this line differs from the reference strain, or generate a genome sequence for this fly strain, 
but we'll get to that later on in the course.) We are also going to use two different (but popular) mapping tools, bwa and bowtie. 
Among their differences is that bowtie (while smokin' fast) does not deal with "gapped" alignments, i.e. it 
does not handle insertion/deletions well. 

Getting the data
----------------

If you just finished the FastQC tutorial, you can keep working on the same machine. Otherwise, launch an EC2 instance, make a volume out of the 
same snapshot that we did before, and mount it. The data for this tutorial are in /data/drosophila::

  %% cd /data/drosophila/
  
For this example, you'll also need the Drosophila reference genome::

  %% curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.37_FB2011_05/fasta/dmel-all-chromosome-r5.37.fasta.gz
  %% gunzip dmel-all-chromosome-r5.37.fasta.gz

Installing and running bwa
--------------------------
  
To actually do the mapping, we need to download and install bwa.

First we are going to grab the source files for bwa from sourceforge, using curl. It is important to 
know that we need to specify a few flags to let the program know that we want to leep the name (and filetype) for 
the file (-O) and that curl should follow relative hyperlinks (-L), to deal with redirection of the file site 
(or else curl won't work with sourceforge).::

  %% curl -O -L http://sourceforge.net/projects/bio-bwa/files/bwa-0.5.9.tar.bz2

Now we want to uncompress the tarball file using "tar". x extracts, v is verbose (telling you what it is doing), 
f skips prompting for each individual file, and j tells it to unzip .bz2 files.::

  %% tar xvfj bwa-0.5.9.tar.bz2
  %% cd bwa-0.5.9

The make command calls a program that helps to automate the compiling process for the program.::

  %% make

(Note: If your system doesn't have make installed already, you'll need to run "apt-get install make" before you 
can build bwa -- but if you are using the right AMI, it should be pre-installed.)

Copy the executable for bwa to a directory for binaries which is in your shell search path::  

  %% cp bwa /usr/local/bin
  %% cd ..

If you want to see what is in your shell search path (which can be modified later), run::

  %% echo $PATH

Now there are several steps involved in mapping our sequence reads and getting the output into a usable form. 
First we need to tell bwa to make an index of the reference genome; this will take a few minutes, so we've already
got the index already generated in the data directory, but if you want to try it yourself, you can run (but see the note below first!!)::

  %% bwa index dmel-all-chromosome-r5.37.fasta

Note: This step takes several minutes. If you run it, it will overwrite the index files we have already made, 
so don't run the above line exactly; instead, create a copy of the reference genome and then index the copy instead, 
so that we can preserve our pre-computed reference index::

  %% cp dmel-all-chromosome-r5.37.fasta dmel-all-chromosome-r5.37.copy.fasta
  %% bwa index dmel-all-chromosome-r5.37.copy.fasta
 
Next, we do the actual mapping. These were paired-end reads, which means that for each DNA fragment, we have sequence 
data from both ends. The sequences are therefore stored in two separate files (one for the data from each end), so we 
have two mapping steps to perform. For now, we'll use bwa's default settings. The files you'll be running this on are 
datasets that have been trimmed down to just the first 1 million sequence reads to speed things up, but at the end you'll 
be able to work with the final product from an analysis of the full dataset that we ran earlier (some of these steps take 
upwards of an hour on the full dataset, but just a couple minutes on the trimmed dataset). Run::

  %% bwa aln dmel-all-chromosome-r5.37.fasta RAL357_1.fastq > RAL357_1.sai
  %% bwa aln dmel-all-chromosome-r5.37.fasta RAL357_2.fastq > RAL357_2.sai

These .sai files aren't very useful to us, so we need to convert them into SAM files. In this step, bwa takes the 
information from the two separate ends of each sequence and combines everything together. Here's how you do it 
(this may take around 10 minutes)::

  %% bwa sampe dmel-all-chromosome-r5.37.fasta RAL357_1.sai RAL357_2.sai RAL357_1.fastq RAL357_2.fastq > RAL357_bwa.sam
  
The SAM file is technically human-readable; take a look at it with::

  %% more RAL357_bwa.sam
  
It's not very easy to understand (if you are really curious about the SAM format, there is a 12-page 
manual at http://samtools.sourceforge.net/SAM1.pdf). For now we'll use bowtie to map the same reads, 
and we'll use another tool to visualize these mappings in a more intuitive way. 

bwa options
-----------

There are several options you can configure in bwa. Probably one of the most important is how many mismatches you 
will allow between a read and a potential mapping location for that location to be considered a match. 
The default is 4% of the read length, but you can set this to be either another proportion of the read length, or a fixed integer. 
For example, if you ran::

  %% bwa aln -n 4 dmel-all-chromosome-r5.37.fasta RAL357_1.fastq > RAL357_1.sai
  
This would do almost the same thing as above, except this time, all locations in the reference genome that contain 
four or fewer mismatches to a given sequence read would be considered a match to that read.

Alternatively, you could do::

  %% bwa aln -n 0.01 dmel-all-chromosome-r5.37.fasta RAL357_1.fastq > RAL357_1.sai
  
This would only allow reads to be mapped to locations at which the reference genome differs by 1% or less from a given read fragment.

If you want to speed things up, you can tell it to run the alignment on multiple threads at once (this will only work 
if your computer has a multi-core processor, which our Amazon image does). To do so, use the -t option to specify the 
number of threads. For example, the following line would run in two simultaneous threads::

  %% bwa aln -t 2 dmel-all-chromosome-r5.37.fasta RAL357_1.fastq > RAL357_1.sai

bwa can also handle single-end reads. The only difference is that you would use samse instead of sampe to generate your SAM file::

  %% bwa samse dmel-all-chromosome-r5.37.fasta RAL357_1.sai RAL357_1.fastq > RAL357_1.sam
  
   
Now let us align our reads using bowtie 
---------------------------------------
(Note: For simplicity we are going to put all of the bowtie related files into the same directory. 
For your own work, you may want to organize your file structure better than we have).

Let's get bowtie from Sourceforge::

  %% curl -O -L http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip

unzip the file, and create a directory for bowtie. In this case, the program is precompiled so it comes as a binary executable::

  %% unzip bowtie-0.12.7-linux-x86_64.zip
  
Change directory::

  %% cd bowtie-0.12.7  

Copy the bowtie files to a directory in you shell search path, and then move back to the parent directory (/data/drosophila)::

  %% cp bowtie bowtie-build bowtie-inspect /usr/local/bin

Let's create a new directory, "drosophila_bowtie" where we are going to place all the bowtie results::

  %% cd ..
  %% mkdir drosophila_bowtie
  %% cd drosophila_bowtie
  
Now we are going to build an index of the Drosophila genome using bowtie just like we did with bwa. The original Drosophila reference genome is in the same location as we used before. Again, we have already performed the indexing step (it takes about 7 minutes), so if you want to try it yourself, index a copy so you don't over-write the one we've pre-run for you::

%%  bowtie-build /data/drosophila/dmel-all-chromosome-r5.37.fasta  drosophila_bowtie  

Now we get to map! We are going to use the default options for bowtie for the moment.  Let's go through this. there are a couple of flags that we have set, since We have paired end reads for these samples, and multiple processors. The general format for bowtie is (don't run this)::

  %% bowtie indexFile fastqFile outputFile

However we have some more details we want to include, so there are a couple of flags that we have to set.
-S means that we want the output in SAM format.
-p 2 is for multithreading (using more than one processor). In this case we have two to use.
-1 -2 tells bowtie that these are paired end reads (the .fastq), and specifies which one is which.
  
This should take 35-40 minutes to run on the full dataset so we'll run it on a trimmed version (should take about 3 minutes; later we'll give you pre-computed results for the full set.)::

  %% bowtie -S -p 2 drosophila_bowtie -1 /data/drosophila/RAL357_1.fastq -2 /data/drosophila/RAL357_2.fastq RAL357_bowtie.sam

You may see warning messages like::

  Warning: Exhausted best-first chunk memory for read SRR018286.1830915 USI-EAS034_2_PE_FC304DDAAXX:8:21:450:1640 length=45/1 (patid 1830914); skipping read

We will talk about some options you can set to deal with this.

The bowtie manual can be found here: http://bowtie-bio.sourceforge.net/manual.shtml

Some additional useful arguments/options (at least for me)
-m  # Suppresses all alignments for a particular read if more than m reportable alignments exist.
-v  # no more than v mismatches in the entire length of the read
-n -l # max number of mismatches in the high quality "seed", which is the the first l base pairs of a read.
-chunkmbs  # number of mb of memory a thread is given to store path. Useful when you get warnings like above
--best # make Bowtie "guarantee" that reported singleton alignments are "best" given the options
--tryhard  # try  hard to find valid alignments, when they exit. VERY SLOW.
 
Processing the output for use with Samtools
-------------------------------------------
  
Even the SAM file isn't very useful unless we can get it into a program that generates more readable output or lets us visualize things in a more intuitive way. For now, we'll get the output into a sorted BAM file so we can look at it using Samtools later.

Download and install Samtools::

  %% cd /data/drosophila
  %% curl -O -L http://sourceforge.net/projects/samtools/files/samtools/0.1.16/samtools-0.1.16.tar.bz2
  %% tar xvfj samtools-0.1.16.tar.bz2
  %% cd samtools-0.1.16
  %% make
  %% cp samtools /usr/local/bin
  %% cd misc/
  %% cp *.pl maq2sam-long maq2sam-short md5fa md5sum-lite wgsim /usr/local/bin/
  %% cd /data/drosophila
  
Like bwa, Samtools also requires us to go through several steps before we have our data in usable form. First, we need to have Samtools generate its own index of the reference genome::

  %% samtools faidx dmel-all-chromosome-r5.37.fasta

Next, we need to convert the SAM file into a BAM file. (A BAM file is just a binary version of a SAM file.)::

  %% samtools import dmel-all-chromosome-r5.37.fasta.fai RAL357_bwa.sam RAL357_bwa.bam

Now, we need to sort the BAM file::

  %% samtools sort RAL357_bwa.bam RAL357_bwa.sorted
  
And last, we need Samtools to index the BAM file::

  %% samtools index RAL357_bwa.sorted.bam

Let us do this again for the bowtie output. We move back to the drosophila_bowtie directory (you could do this all from the other directory, but it gets harder to read the command with long pathnames)::

  %% cd drosophila_bowtie
  %% samtools import ../dmel-all-chromosome-r5.37.fasta.fai RAL357_bowtie.sam RAL357_bowtie.bam

Now, we need to sort the BAM file (also slow)::

  %% samtools sort RAL357_bowtie.bam RAL357_bowtie.sorted
  
And last, we need Samtools to index the BAM file::

  %% samtools index RAL357_bowtie.sorted.bam
  
All done! Now we can use the sorted BAM file in Samtools to visualize our mappings, generate lists of SNPs, and call consensus sequences. We'll get to all of that later on today and in the rest of the course.

Viewing the output with TView
-----------------------------

Before we can use TView to compare Bowtie and BWA mappings, we need to sort the Bowtie BAM file, and generate an index for it.::

  %% cd drosophila_bowtie
  %% samtools sort RAL357_full_bowtie.bam RAL357_full_bowtie.sorted
  %% samtools index RAL357_full_bowtie.sorted.bam

Now that we've generated the files, we can view the output with TView. We'll compare two different sorted::

  %% cd ..
  %% samtools tview RAL357_full_bwa.sorted.bam

Now open an additional terminal window, and load the Bowtie mapping file there as well.::
  
  %% cd /data/drosophila/drosophila_bowtie
  %% samtools tview RAL357_full_bowtie.sorted.bam ../dmel-all-chromosome-r5.37.fasta

To view the tview help, type '?'