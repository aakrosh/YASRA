#for a full scale assisted assembly:
#  keep a FASTA/FASTQ file with the reads from the target genome,
#  keep a FASTA file with the reference genome sequence,
#  change the values of the variable C,READS,REFERENCE and then type:
#
#  make TYPE= ORIENT= PID= &> Summary
#  
#  The choices for TYPE are '454' and 'solexa'. '454' refers to reads from 454
#  technology which are more than 100 bp long. 'solexa' refers to shorter reads
#  in base-space from Illumina or SOLiD.
#	
#  The choices for ORIENT are 'linear' and 'circular'. It refers to the 
#  orientation of the reference genome. (example: mtDNA would be circular)
#
# The choices for PID are:
# 
# --for 454 reads:
#    'same'       : about 98% identity
#    'high'       : about 95% identity
#    'medium      : about 90% identity
#    'low'        : about 85% identity
#    'verylow'    : about 75% identity
#    'desperate'  : realy low identity (rather slow)
#
# --for solexa reads:
#	'same'        : about 95% identity
#	'medium'      : about 85% identity
#	'desperate'   : low scores (rather slow)
#
# The module "best_hit" selects one alignment for each read. The user can select
# one of the following options:
#
# -u : Ignore reads with multiple alignments. This is the default behavior in
#      this Makefile.
# -s : Choose the place with the highest number of matches. If two alignments 
#      have equal number of matches, then choose one of them randomly.
# -b : Choose the best alignment only if it has x% (x can be user-specified, e.g
#      -b 3 for x=3) more matches than the second best hit. (We use x=3 and use
#      this option internally for whole genome analyses).
#
# The recursive  mode of assembly is rather slow because we try to walk through
# the small gaps and close the contigs. However, if we want to ignore the gaps 
# and want quicker results (and this is the mode I would advise most of the 
# time) is by typing:
#
#  make single_step TYPE= ORIENT= PID= &> Summary
#
# Final output of the process are the following files:
#   Final_Assembly : a FASTA file with the consensus sequence
#   alignments.sam : a SAM file with the alignments. This can be converted to a
#                    BAM file and view with a multitude of viewers
#   contigs.ace    : a  ACE file with the alignments for viewing with Hawkeye
#
# The file alignments.sam replaces "Assembly.qual" in the previous versions of
# YASRA. alignments.sam can be used to get all the information present in
# Assembly.qual by use of SAMtools (http://samtools.sourceforge.net/).
#
# For more information regarding options and tools in YASRA please contact :
# 
# ratan@bx.psu.edu
#
C=./../src

# the fasta file with the reads
READS=454.fa

# the fasta file with the reference sequence
TEMPLATE=rhino_template.fa

# the orientation of the reference sequence
ORIENT=linear

# the type of input sequence 454/solexa
TYPE=454

# the expected percent identity between the reference and target genomes
PID=same

ifeq ($(TYPE), 454)
	MAKE_TEMPLATE=min=150
    ifeq ($(PID),same)
            Q=--yasra98
    endif
    ifeq ($(PID),high)
            Q=--yasra95
    endif
    ifeq ($(PID),medium)
            Q=--yasra90
    endif
    ifeq ($(PID), low)
            Q=--yasra85
    endif
    ifeq ($(PID), verylow)
            Q=--yasra75
    endif
    ifeq ($(PID), desperate)
            Q=Y=2000 K=2200 L=3000
    endif
	R=--yasra98
	names=darkspace
endif

ifeq ($(TYPE), solexa)
	MAKE_TEMPLATE=N=100 min=30
	SOLEXA=-solexa
    ifeq ($(PID),same)
            Q=--yasra95short
    endif
    ifeq ($(PID),medium)
            Q=--yasra85short
    endif
    ifeq ($(PID), desperate)
            Q=T=0 W=7 K=1200 L=1400
    endif
	R=--yasra95short
	names=full
endif

ifeq ($(ORIENT), circular)
	CIRCULAR=--circular
endif

all:step1 step2 step3 step4 step5

step1:
	# Assemble on the original template
	make assemble_hits T=$(TEMPLATE) P="$Q" V=70 N=1

step2:
	$C/genomewalker Assembly1 hits_$(TEMPLATE) $(TEMPLATE)
	@rm Assembly* hits* template[0-9]* 

step3:
	# find if there are reads that align without the coverage option
	lastz template[multi] $(READS)[nameparse=${names}] $R \
		--coverage=70 --ambiguous=iupac \
		--format=general:name1,zstart1,end1,text1,name2,strand2,zstart2,end2,text2,nucs2 |\
	$C/best_hit -u |\
	sort -k 1,1 -k 2,2n -k 3,3n > hits_template
	lastz template[multi] $(READS)[nameparse=${names}] --yasra85short \
		--coverage=50 --ambiguous=iupac \
		--format=general:name1,zstart1,end1,text1,name2,strand2,zstart2,end2,text2,nucs2 |\
	$C/best_hit -u |\
	sort -k 1,1 -k 2,2n -k 3,3n > rejects_template

step4:
	$C/trim_assembly template hits_template rejects_template > AssemblyX
	$C/make_template AssemblyX noends $(MAKE_TEMPLATE) > final_template
	@rm AssemblyX template 
	@rm hits_template rejects_template
	
step5:
	 make final_assembly T=final_template P="$R" V=70
	@rm final_template 
	
stepx:
	$C/make_template Assembly$W $(MAKE_TEMPLATE)> template$X
	make assemble_hits T=template$X P="$R" V=70 N=$X

assemble_hits:
	lastz $T[multi] $(READS)[nameparse=${names}] $P \
		--coverage=$V --ambiguous=iupac \
		--format=general:name1,zstart1,end1,text1,name2,strand2,zstart2,end2,text2,nucs2 |\
	$C/best_hit -u |\
	sort -k 1,1 -k 2,2n -k 3,3n > hits_$T
	$C/assembler -o -c -h hits_$T > Assembly_$N
	$C/genomewelder $(CIRCULAR) -n 10 -i 95 -c 95 Assembly_$N > Assembly$N

single_step:
	make final_assembly T=$(TEMPLATE) P="$Q" V=70

final_assembly:
	lastz $T[multi] $(READS)[nameparse=${names}] $P \
		--coverage=$V --ambiguous=iupac \
		--format=general:name1,zstart1,end1,text1,name2,strand2,zstart2,end2,text2,nucs2 |\
	$C/best_hit -u |\
	sort -k 1,1 -k 2,2n -k 3,3n |\
	$C/assembler -r -o -c \
		-h /dev/stdin \
		-s alignments.sam \
		-a contigs.ace > Final_Assembly

clean:
	rm Final_Assembly alignments.sam contigs.ace
