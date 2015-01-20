YASRA
=====

## REQUIREMENTS
YASRA should work on any standard 64 bit Linux environment with gcc and python.
However we have tested the compilation on the following platforms (output of
gcc --version and uname -a):

a) gcc (GCC) 4.1.2 20080704 (Red Hat 4.1.2-50)
   Linux 2.6.18-164.el5 #1 SMP x86_64 GNU/Linux

YASRA uses LASTZ (http://bx.psu.edu/miller_lab for released version and
http://www.bx.psu.edu/~rsharris/lastz/newer for newer version) for aligning the
sequences to the reference genome. Please install LASTZ (the newest version on
http://www.bx.psu.edu/~rsharris/lastz/newer) and add the LASTZ binary in your
executable/binary search path before installing YASRA. 

## SUMMARY
YASRA (Yet Another Short Read Assembler) performs comparative assembly of short
reads using a reference genome, which can differ substantially from the genome
being sequenced. Mapping reads to reference genomes makes use of LASTZ (Harris
et al), a pairwise sequence aligner compatible with BLASTZ. Special
scoring sets were derived to improve the performance, both in runtime and
quality for 454 and Illumina sequence reads. 

## INSTALLATION
The programs can be compiled by following the following recipe:

  % configure --prefix=/usr/local (or whatever installation path you prefer)
  % make
  % make install

This complies all the components of the pipeline and puts the binaries in the
folder $prefix/bin.  For more in depth instructions, consult the INSTALL file.

Please add the $prefix/bin folder to your executable/binary search path to
complete the installation.

## DESCRIPTION
Our proposed pipeline proceeds in several steps. The first step involves 
aligning the reads to a distant template followed by trimming of the reads at 
the ends. Subsequently the reads are assembled to form the consensus sequence.
This process can be repeated to close some of the gaps in the assembled
sequence. The last iteration is followed by an error-detection module, which
aims to detect areas with low coverage or with possible assembly errors. 

## READ NAMES
Up till YASRA version 2.31 only the first space delimited token (in the read
header) was considered as the unique read identifier. For example if a 454 read
had the header:

>GAOCPZU01A7HGZ length=115 xy=0378_0673 region=1 run=R_2010_01_14_15_16_37_

I considered GAOCPZU01A7HGZ as the unique read identifier. 

Similarly, for an Illumina read downloaded from the Short Read Archive with the
header:

@HWI-EAS313:7:1:8:1009#0/1 length=76

HWI-EAS313:7:1:8:1009#0/1 was considered to be the unique read identifier. 
But recent Illumina runs use read headers such as 

@HWI-ST978:170:D1032ACXX:5:1101:1461:2178 1:N:0:CGGAAT

where the first token is not unique. In order to handle such cases all tokens in
Illumina reads are now considered in the read names from YASRA version 2.32
onwards.  This behavor is turned on in the Makefile in the line 116 of the test
dataset:

names=full

For the 454 datasets, I still consider only the first token as the read
identifier. The same behavior (as Illumina, where all tokens are considered) can
be  turned on for 454  sequences as well, by changing line 100 in the Makefile
from

names=darkspace

to

names=full
 

## TEST-DATASET
A sample toy dataset can be found in the "test_data" sub-directory. The reads
(454.fa) are from a Black rhinoceros sample. We use the non VNTR region of the Indian rhinoceros to construct the consensus sequence of the mtDNA of the Black
rhinoceros. Please take a look at the Makefile in the directory and use it to
assemble the mtDNA. 

The Makefile can be modified to suit your needs for a particular project. In
summary YASRA has two modes:
a) single_step : where the reads are aligned to the reference sequence and a
   consensus sequence is constructed out of that. This is the preffered method
   that should be used in most cases.
b) recursive : where the reads are aligned to the reference and a consensus
   sequence is generated. That sequence is then used as a reference and this is
   continued till we cannot improve the assembly any further. Please use this
   mode with caution. I would always recommend using the single_step mode before
   you play around with the recursive mode for YASRA.

For every project, create a copy of the Makefile in the test_data directory and
modify the following variables in the Makefile:

C        : this should point to the directory with all the binaries for YASRA.
READS    : this should point to the fasta file with the reads for the project.
TEMPLATE : this should point to the fasta file with the reference genome.
ORIENT   : linear/circular depending on whether the reference genome is linear 
           or circular. For example mtDNA is circular.
TYPE     : 454/solexa depending on whether the reads were sequenced using 454 
           or Solexa/Illumina.
PID      : for 454 reads this could be
                'same'       : about 98% identity between target & reference
                'high'       : about 95% identity between target & reference
                'medium      : about 90% identity between target & reference
                'low'        : about 85% identity between target & reference
                'verylow'    : about 75% identity between target & reference
                'desperate'  : realy low identity (rather slow)
           for Solexa/Illumina reads this could be 
                'same'        : about 95% identity between target & reference
                'medium'      : about 85% identity between target & reference
                'desperate'   : low scores (rather slow)
 
The user can override these options from the command-line. For example if the
Makefile had PID=same, and an user wanted to attempt to run the pipeline with
PID=low, then they could modify the Makefile or just type on the command line:

    make PID=low

Please ensure that you have the latest version of LASTZ installed. You can check
the version installed by using the following command : 

    lastz --version


Some of the parameters can only be changed in the Makefile. The defaults should
work for most cases, but users might want to make these changes for their
purpose for specific projects.

These are the options for the module "best_hit" which selects one alignment per
qualifying read. The user should use only one of the following options:

-u : Ignore reads with multiple alignments. This is the default behavior in
     the Makefile with the test dataset.
-s : Choose the place with the highest number of matches. If two alignments 
     have equal number of matches, then choose one of them randomly.
-b : Choose the best alignment only if it has x% (x is user-specified, e.g
     -b 3 for x=3) more matches than the second best hit. (We use x=3 
     internally for whole genome analyses).
 
In cases where the the assembly does not benefit from the iterative process (if
that option was chosen), it exists with the message

```
This assembly does not benefit from the recursive process. Please run "make
single_step" instead.
make: *** [step2] Error 1
```

Please run 

make clean 

and then run the single_step assembly.

## IMPROVING THE ASSEMBLY

After the initial assembly is done, one of the ways to visualize it is using LAJ
(http://globin.bx.psu.edu/dist/laj). If your reference sequence is reference.fa
(which is a single sequence) and the final assembly is Final_Assembly (can be multiple contigs), then do the following:

lastz reference.fa Final_Assembly --chain > fake.bz

laj fake.bz

This should bring up a dot plot of the alignments between the reference and the
assembled contigs. One can identify the regions that can be improved from these
alignments. Specifically, I look at these alignments to identify neighboring
contigs that can be merged together. 


## DETAILS
The most important component of the pipeline is a binary called "assembler". It
requires the following arguments:

-d, --debug : If set, print out debug information along with the results. This 
              is much slower than the actual version and should only be run 
              for debugging. [Default : not set]
-h, --hits  : This is the only required argument for the module. The file should
              be the result of running LASTZ to align reads to a reference
              sequence with the output              
              format=general:name1,zstart1,end1,name2,strand2,zstart2,end2,nucs2
-a, --ace   : Name of the output ACE file for the assembly [Default : NULL]
-s, --sam   : Name of the output SAM file for the assembly [Default : NULL]
-m, --cov   : Only print out contigs where the average depth of coverage is
              greater than or equal to this value [Default : 0 ]

-o, --orient   : Does the sequence in the LASTZ alignment need to be oriented
                 prior to assembly?
-c, --declone  : Throw away putative PCR duplicates when assembling the contigs.
                 A putative PCR duplicate is a read that has the same start 
                 and end points on the reference as another read.
-r, --realign  : Realign the reads after the initial consensus is made, to
                 improve the alignments and hence the consensus. This uses the 
                 algorithm described in Anson et al. The iterations continue if
                 an objective score can be improved in the subsequent iteration.
-1, --realign1 : Just run the realignment for 1 iteration.
