Title: Installing and Running NCBI BLAST
Tags: blast tutorial MSU-NGS-2011
##rest

You should start this tutorial at a prompt that looks something like this::

   root@ip-10-82-233-6:~#

Type 'cd' to go to your home directory on your EC2 machine.

Now, use your Web browser on your laptop to go to:

   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST

right- or control-click on the file ending with 'x64-linux.tar.gz', and
"copy link URL".  This is the file for 64-bit (large) Linux machines, which
is what our EC2 instance is.  (The current URL is: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.25+-x64-linux.tar.gz)

Now use the 'curl' program to download it to your Amazon computer::

 %% curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.2.25+-x64-linux.tar.gz

Here, 'curl' is a program that takes a Web link and downloads it via the
command line; in this case, it's grabbing that file and saving it into your
current directory.

After it completes, you should see the file in your local directory::

 %% ls ncbi-*.tar.gz

This is a .tar.gz file, which is kind of like a zip file.  You need to
use the 'tar' program to unpack it (you could use 'unzip' if it were a .zip
file)::

 %% tar xzf ncbi-*.tar.gz

This will create a new subdirectory, 'ncbi-blast-2.2.25+'::

 %% ls
 Dropbox  ncbi-blast-2.2.25+  ncbi-blast-2.2.25+-x64-linux.tar.gz

If you look in the blast subdirectory, you will see a few more files, most of
which are directories::

 %% ls ncbi-blast-2.2.25+
 bin  ChangeLog  doc  LICENSE  ncbi_package_info  README

In this case, we want to put everything in that bin/ directory into
a common place where UNIX knows to look for programs to run.  One such
place (that, by convention, is a good place to install things that don't
come with the computer) is /usr/local/bin::

 %% cp ncbi-blast-2.2.25+/bin/* /usr/local/bin

Now, let's go to a new section of the machine. ::

 %% cd /mnt

This goes to the folder named '/mnt', which is on another (bigger)
disk.  We'll explain this more tomorrow.

Now lets grab some biggish files to work with... the mouse and
zebrafish reference proteomes!

Go to ftp://ftp.ncbi.nlm.nih.gov/refseq/ in your browser and explore a
bit.  You'll see there's a bunch of files and directories; in this
case, we want to go grab the mouse and zebrafish protein sets. So,
grab
ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz
and
ftp://ftp.ncbi.nlm.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.protein.faa.gz::

 %% curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz
 %% curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.protein.faa.gz

These files aren't .tar.gz files, they're just .gz files -- the .faa means
"fasta".  'gz' is a compression scheme for single files; to get at the
contents, do uncompress both of them with this command: ::

 %% gunzip *.gz

If you use 'ls', you'll see that the files have turned into 'mouse.protein.faa'
and 'zebrafish.protein.faa'::

 %% ls

You can also take a look at the contents of the files with the 'more'
program, which pages through the files. ::

 %% more mouse.protein.faa

Use the spacebar to scroll down, and 'q' to exit before reaching the
end of the file.  You can also look at the zebrafish file::

 %% more zebrafish.protein.faa

Now, let's convert them into BLAST databases::

 %% makeblastdb -in mouse.protein.faa -dbtype prot
 %% makeblastdb -in zebrafish.protein.faa -dbtype prot

This lets us use BLAST to query the databases for matches.

Before we do a *big* BLAST, let's start by doing a small one, just to
check that it's all working.  To do that, we'll skim off some
sequences from the top of the file::

 %% head zebrafish.protein.faa

The problem here is that 'head' by default only selects the first 10
lines of a file, which may not be a complete set of FASTA records -- so you
may have to tweak things.  In this case, the first 14 lines are complete::

 %% head -14 zebrafish.protein.faa

Let's take the output of 'head' and put it in a file, 'zebrafish.top',
that we can use for other purposes::

 %% head -14 zebrafish.protein.faa > zebrafish.top

OK, great!  Now let's run a BLASTP comparing these zebrafish sequences
to the mouse proteins, and we'll put the results in a file 'xxx.txt'::

 %% blastp -query zebrafish.top -db mouse.protein.faa -out xxx.txt

(The file name 'xxx.txt' is just a throwaway file name, something
that we can look at and see is a test file.  You can use your own
convention; I usually go with something short and recognizably
silly, like 'xxx', 'yyy', 'foo', etc.)

OK, now take a look at that file with 'more'::

 %% more xxx.txt

Yep, looks like BLAST output to me!

There's all sorts of things you can do to alter the BLAST output; run
'blastp' to get a list of those options.  For example, '-evalue 1e-6'
will set the e-value cutoff at 1e-6, above which nothing will be
displayed.

Now let's run a bigger BLAST, all zebrafish proteins against all mouse
proteins::

 %% blastp -query zebrafish.protein.faa -db mouse.protein.faa -out zebrafish.x.mouse &

This is going to take a while, which is why we told the computer to
give us back a command prompt while blastp runs (that's what the &
does).

So, how long is it going to take?  We can guesstimate by looking at
how many sequences have been processed since we started.  To do that, run ::

 %% grep Query= zebrafish.x.mouse

OK, that gives us all the query lines -- now what?  Let's count them with
'wc -l'::

 %% grep Query= zebrafish.x.mouse | wc -l

Here, | is what's known as a 'pipe', telling the command line to take
the output of 'grep' and send it to the command 'wc', which counts
words, lines, and paragraphs.  The '-l' tells wc to count the lines
only.

Compare that number to the number of sequences in the zebrafish protein database::

 %% grep ^'>' zebrafish.protein.faa | more

to see the FASTA headers, and ::

 %% grep ^'>' zebrafish.protein.faa | wc -l

to count all the sequences.

Last, but not least -- let's run a quick script to convert the file into
a set of CSV matches::

 %% python ~/Dropbox/ngs-scripts/blast/blast-to-csv.py zebrafish.x.mouse > ~/Dropbox/zebrafish-mouse.csv

Take a look at the script and see if you can understand what it does::

 %% more ~/Dropbox/ngs-scripts/blast/blast-to-csv.py

Before you leave for lunch:
---------------------------

Let's start a *second* BLAST, all of mouse against all of zebrafish::

  %% blastp -query mouse.protein.faa -db zebrafish.protein.faa -out mouse.x.zebrafish &

...now the computer can work while we eat!

When we come back, we can work through a reciprocal BLAST example.

.. @@ save these files

.. # scripting etc.
.. @@ install biopython, ez_seutp, blastkit??