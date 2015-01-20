YASRA
======
Reference based assembler

## REQUIREMENTS
YASRA uses LASTZ to map the reads on to a reference. LASTZ can be freely
downloaded from http://www.bx.psu.edu/miller_lab/ YASRA should work on any
standard 64 bit Linux environment with gcc, python and gnuplot.

## INSTALLATION
Please install LASTZ and ensure that it is included in your PATH. All the
remaining modules for YASRA are included in this directory. Simply follow the
instructions in the INSTALL file and you should be able to install the binaries
in a folder of your choice

## NOTES
### Comparative Assembler Pipeline:
The requirements of this module are the following:
1) Given a set of hits(the reads and the positions of the hit on the template),
this assembler must be able to assemble those reads to form the contigs.
2) The assembler output must be displayable in the standard ACE format as well
as a format acceptable to Realigner, since that will be a part of the pipeline. 
3) The assembler should be able to handle quality values, even though there
might not be much that can be done with them as of now. Hence, the assembler
must be able to read the quality value files, and the quality values must be a
part of the data structure.
4) The assembler should work for both 454 and Solexa reads and the decision
whether it is solexa or 454 should be made with a command line switch on
runtime.
 
### ACE Format
The .ACE format is produced by phrap as well as by most other assemblers (including Arachne, TIGR Assembler, CAP, etc.)

Example:

CO 1 30502 510 273 U
CCTCTCC*GTAGAGTTCAACCGAAGCCGGTAGAGTTTTATCACCCCTCCC

BQ
 20 20 20 20 20 20 20 20 20 20 20 20 20

 AF TBEOG48.y1 C 1

 BS 1 137 TBEOG48.y1

 RD TBEOG48.y1 619 0 0
 CCTCTCC*GTAGAGTTCAACCGAAGCCGGTAGAGTTTTATCACCCCTCCC

 QA 1 619 1 619

Contig identifiers (starting with CO) list the IDs  (1 in the example),
the number of bases (30502), number of reads (510), and number of "base
segments" (273) as well as whether the contigs is in the forward orientation
(Uncomplemented) or reversed (Complemented).  In general, the output of an
assembler has all contigs listed as "U".

The consensus sequence is padded, the gaps being represented as *s
instead of dashes, and follows immediately after the CO line.
Consensus quality values are provided for the bases alone, the gaps not 
being represented. These quality values, in phred-like format follow
immediately after the BQ line. 

The AF lines (one per aligned read) contain information of
whether the read is complemented (C) or not (U) followed by a 1-based offset
in the consensus sequence.  Note that the offset refers to the beginning of
the entire read in the alignment, not just the clear range.  Thus the read
acaggATTGA will have an offset of 1 even though the consensus truly starts at
position 6.

The BS lines indicate which read was used to calculate the consensus 
between the specified coordinates.  These lines can, in general, be ignored 
as they are an artifact of the algorithms used to compute the consensus 
sequence.

The sequence of each read is explicitly provided after each RD line.  The 
sequence is padded with *s and is already complemented if necessary.
The QA line following each read contains two 1-based ranges. The second 
range represents the clear range of the read, with respect to the read 
sequence (padded and potentially complemented) as provided in the RD record.

### Design of the data structures

The basic design of the stuff used in the assembler is in contig.h.
A 'contig' refers to an assembled contig by the assembler. The 'next' member of
the contig refers to the next contig in the assembly list.

Contig: is formed of columns, which are in a list.

------------------->

.-----------.-----------.-----------.----------.----------.----------.
|    A      |     C     |     C     |    G     |    G     |     A    | Read 1
|    A      |     C     |     C     |    G     |    G     |----------' Read 2
|    A      |     C     |     C     |    G     |    G     |	           Read 3
|    A      |     C     |     C     |    G     |    G     |            Read 4
|    A      |     C     |     C     |    G     |    G     |              |
|    A      |     C     |     C     |    G     |    G     |              |
|    A      |     C     |     C     |    G     |    G     |              |
|    A      |     C     |     C     |    G     |    G     |
|    A      |     C     |     C     |    G     |    G     |
|    A      |     C     |     C     |    G     |    G     |
|    C      |     C     |     C     '----------'----------'
|    C      |     C     |     C     |
|    C      |     G     |     C     |
|    C      |     C     |     C     |
'-----------'-----------'-----------'

Each base in a column is called an element.

### Deletions in the target compared to the reference
If the deletion in the target is less than the length of a read, it is taken
care of by the assembler. The trimmer assigns it a long mapped length on the
reference, which makes it a non-contained read, leading it to bridge the gap.
However if the deletion in the target is higher than the length of a read then
this cannot be bridged in a single iteration. Hence it  would lead to a break in
the contigs.

### Insertion in the target compared to the reference
Similarly if the insertion in  the target genome is smaller than the length of a
read, then the assembler is able to bridge that area. However if the length of
the read is greater than a read, then again, this cannot be bridged with just
one iteration.

### LAV format
The details of the lav format can be found in the LASTZ documentation.

## TEST-DATASET
Some additional tools are provided with the pipeline in tools.mk. They are
simple awk commands that can be helpful in looking at different aspects of the
assembly. 

A sample toy dataset can be found in the "test_data" sub-directory. The reads
are a sample of 10X coverage from the sequence correct.fa. They were generated
using ReadSim (http://www-ab.informatik.uni-tuebingen.de/software/readsim/) with
a mean length of 100 bp. The directory contains a Makefile, which demonstrates
the use of the pipeline in steps, for a reference which is approximately 95.6%
similar to the correct sequence. That Makefile discusses one iteration of the
assembler on the reference to generate contigs. It also illustrates the
conversion of the ACE file into a bank which can then be viewed using Hawkeye.

A complete Makefile which uses the iterative mode of the pipeline can be made 
as shown at the end of this document. 

```
# Put the reads file as 454.fa, the reference as reference.fa and run this
# makefile as
#	make TYPE=454 ORIENT=linear PID=same
#	
# Options for TYPE are 454,solexa
# Options for ORIENT are linear,circular
# Options for PID are same,high,medium,low,verylow and denote percent identity
#     to the reference. Same 98%, High 95% and so on.

#this needs to be changed to point to the directory with the binaries.
C=

#this is the length of the ids of the reads. 
WL=`cat 454.fa | grep '>' | awk '{print length($$1)}' | sort -nr | head -1`

#is this 454 or solexa data
TYPE=454

#is this a circular or a linear genome
ORIENT=linear

#maximum length of a read
MAX=`cat 454.fa | awk '{if(substr($$0,1,1) == ">"){if(a>max){max=a}; a=0} else{a = a+length($$0)}}; END{print max}'`

# Q gives parameters for finding weak matches to a rather distant species
# R gives parameters for finding high-identity matches to endogenous DNA

ifeq ($(TYPE), 454)
	MAKE_TEMPLATE=min=150
	ifeq ($(PID),same)
		Q=98
	endif
	ifeq ($(PID),high)
		Q=95
	endif
	ifeq ($(PID),medium)
		Q=90
	endif
	ifeq ($(PID), low)
		Q=85
	endif
	ifeq ($(PID), verylow)
		Q=75
	endif
	R=98 
endif

ifeq ($(TYPE), solexa)
	MAKE_TEMPLATE=N=100 min=30
	SOLEXA_INSERT=N=100
	SOLEXA=-solexa
	ifeq ($(PID),same)
		Q=95short
	endif
	ifeq ($(PID),high)
		Q=95short 
	endif
	ifeq ($(PID),medium)
		Q=85short
	endif
	ifeq ($(PID), low)
		Q=85short
	endif
	ifeq ($(PID), verylow)
		Q=85short
	endif
	R=95short 
endif

ifeq ($(ORIENT), circular)
	CIRCULAR=-circular
endif

TEMPLATE=reference.fa
COMPARE=$(TEMPLATE)

CON_INFO=info.txt
REJECT=hits_reject.fa
REJ=Rejects.txt
REP=repeats.txt
HITS=hits_final.fa

all:step1 step2 step3 step4 step5

single_step:
	make final_assembly T=$(TEMPLATE) V=60 P="$Q" I=80 S="-sig=1"
	rm MAlign ReMAlign

transcriptome:
	$C/addX genes | grep -v '>' > $(TEMPLATE)
	make final_assembly T=$(TEMPLATE) V=60 P="$Q" I=60 S="-sig=1"
	rm MAlign ReMAlign

step1 :
	#Assemble on the original template:
	make assemble_hits T=$(TEMPLATE) V=60 I=70 P="$Q" N=1

step2:
	touch fake.txt 
	$C/finisher $(REJ) $(REP)
	@rm Assembly* hits* template[0-9]* $(REP) fake.txt

step3:
	#Determine difference between reads and assembly:
	lastz template 454.fa --yasra$R | \
	$C/lav_processor | \
	$C/template_hits template 454.fa cov=70 pct=90 stats \
	   discard=$(REJECT) > read_hits
	make plot H=read_hits

step4:
	#Trim Assembly:
	if [ -s plot_rejects.eps ]; then \
		mv plot_rejects.eps old_plot_rejects.eps; \
	fi
	$C/trim_assembly template read_hits $(REJECT) min=1 > AssemblyX
	$C/make_template AssemblyX noends $(MAKE_TEMPLATE) info=$(CON_INFO) \
		> ftemplate
	make assemble_hits T=ftemplate V=70 I=90 P="$R" N=Y
	$C/make_template AssemblyY noends $(MAKE_TEMPLATE) info=$(CON_INFO) \
		> final_template
	lastz final_template 454.fa --yasra$R | \
	$C/lav_processor | \
	$C/template_hits final_template 454.fa cov=70 pct=90 stats \
	   discard=$(REJECT) > read_hits
	make plot H=read_hits
	@rm read_hits AssemblyX AssemblyY ftemplate Assembly_ftemplate \
		hits_ftemplate.fa

step5:
	make final_assembly T=final_template V=70 P="$R" I=90 S="-sig=1"
	@rm final_template template MAlign ReMAlign

final_assembly:
	time lastz $T 454.fa --yasra$P | \
	time $C/lav_processor $S  | \
	time $C/template_hits $T 454.fa cov=$V pct=$I \
		discard=$(REJECT)  > $(HITS)
	cat $(HITS) | grep '>Hit' | awk '{print $$4, $$1}' | sort -n |  \
		uniq -w $(WL) -D | awk '{print $$2}' > $(REP)
	time $C/assembler hits_final.fa -ace=Contigs.ace -rejects=$(REJ) \
		-repeats=$(REP) $(SOLEXA) -max_length=$(MAX) > MAlign
	time $C/realigner -ace=Contigs.ace < MAlign > ReMAlign
	time $C/consense -ace=Contigs.ace -amb=Amb.txt -profile=Assembly.qual \
	 < ReMAlign > Final_Assembly
	toAmos -ace Contigs.ace -o  - | bank-transact  -f -b bank -m -

plot:
	make -f $C/tools.mk rejects H=$H R=$(REJECT) C=$(CON_INFO);

stepx:
	#
	#Assemble on the endogenous contigs:
	cp fake.txt $(CON_INFO)
	$C/make_template Assembly$W info=fake.txt $(MAKE_TEMPLATE)> template$X
	make assemble_hits T=template$X V=70 P="$R" I=90 N=$X


assemble_hits : 
	time lastz $T 454.fa --yasra$P | \
	$C/lav_processor  | \
	$C/template_hits $T 454.fa cov=$V pct=$I discard=$(REJECT) \
		stats > hits_$T.fa
	time make assemble HITS=hits_$T.fa
	$C/welder Assembly_$T $(CIRCULAR) $(SOLEXA) | $C/consense > Assembly$N

assemble:
	cat $(HITS) | grep '>Hit' | awk '{print $$4, $$1}' | sort -n |  \
		uniq -w $(WL) -D | awk '{print $$2}' > $(REP)
	$C/assembler $(HITS) -repeats=$(REP) $(SOLEXA) \
		-rejects=$(REJ) -max_length=$(MAX) | \
	$C/realigner > ReMAlign
	$C/consense -amb=Amb.txt -profile=Assembly.qual \
	    < ReMAlign > Assembly_$T

#optionally fix the assembly for homopolymer runs, and areas where we are not 
#sure. This should only be use if the reference and the assembled sequence
#are really close or the same.
fix:
	$C/fix_asm Final_Assembly subs.txt Amb.txt Final_Assembly -H -S \
		-E > Fixed_Final_Assembly
	mv Fixed_Final_Assembly Final_Assembly
```
