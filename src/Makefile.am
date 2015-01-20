## Makefile.am for YASRA

if DEBUG
OPTIMIZATIONS = -ggdb -DDEBUG
else
if PROFILE
OPTIMIZATIONS = -ggdb -pg -DNDEBUG -O2 -funroll-all-loops $(ARCH)
else
if TESTING
OPTIMIZATIONS = -ggdb
else
OPTIMIZATIONS = -DNDEBUG -O2 -funroll-all-loops $(ARCH)
endif
endif
endif

AM_CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall -Werror \
            -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
            -Wcast-align -Winline -Wnested-externs -Wextra $(OPTIMIZATIONS) \
            -std=c99

###############################################################################

bin_PROGRAMS = addX make_template template_hits lav_processor assembler \
               realigner consense substitutions trim_assembly welder fix_asm

bin_SCRIPTS  = finisher tools.mk

###############################################################################

UTIL_SRC   = util.h util.c
CHR_SRC    = charvec.h charvec.c
SEQ_SRC    = charvec.h charvec.c seq.h seq.c
SLIST_SRC  = slist.h slist.c
HASH_SRC   = hash.h hash.c
CONTIG_SRC = contig.h contig.c
GRAPH_SRC  = graph.h graph.c
IUPAC_SRC  = iupac.h iupac.c

addX_SOURCES = $(UTIL_SRC) $(SEQ_SRC) addX.c

make_template_SOURCES = $(UTIL_SRC) $(SEQ_SRC) make_template.c

template_hits_SOURCES = $(UTIL_SRC) $(SEQ_SRC) template_hits.c

lav_processor_SOURCES = $(UTIL_SRC) $(SLIST_SRC) lav_processor.c

assembler_SOURCES = $(UTIL_SRC) $(SEQ_SRC) $(SLIST_SRC) $(CONTIG_SRC) \
                    $(HASH_SRC) $(GRAPH_SRC) assembler.c

realigner_SOURCES = $(UTIL_SRC) realigner.c

consense_SOURCES  = $(UTIL_SRC) $(CHR_SRC) $(SLIST_SRC) $(IUPAC_SRC) consense.c

substitutions_SOURCES = $(UTIL_SRC) $(SEQ_SRC) substitutions.c

trim_assembly_SOURCES = $(UTIL_SRC) $(SEQ_SRC) trim_assembly.c

welder_SOURCES = $(UTIL_SRC) $(SEQ_SRC) $(SLIST_SRC) $(CONTIG_SRC) welder.c

fix_asm_SOURCES = $(UTIL_SRC) $(SEQ_SRC) $(SLIST_SRC) $(HASH_SRC) fix_asm.c


