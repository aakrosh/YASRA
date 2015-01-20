# --enable-profile
# --enable-debug
# --enable-test

## i386 and x86-64 options: -minline-all-stringops -momit-leaf-frame-pointer

if DEBUG
OPTIMIZATIONS = -ggdb -DDEBUG
else
if PROFILE
OPTIMIZATIONS = -ggdb -pg -DNDEBUG -O2 -funroll-all-loops $(ARCH)
else
OPTIMIZATIONS = -DNDEBUG -O2 -funroll-all-loops $(ARCH)
endif
endif


AM_CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall -Werror \
            -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
            -Wcast-align -Winline -Wnested-externs -Wextra $(OPTIMIZATIONS) \
            -std=c99

################################################################################

bin_PROGRAMS = assembler genomewelder make_template

################################################################################

ASSERTS_SRC = asserts.h
COMMON_SRC  = common.c common.h
CONTIG_SRC  = contig.c contig.h
GRAPH_SRC   = graph.c graph.h
HASHTABLE_SRC = hashtable.c hashtable.h
RUNTIME_SRC = runtime.c runtime.h
SEQUENCES_SRC = sequences.c sequences.h
SLINKLIST_SRC = slinklist.c slinklist.h
UTILITIES_SRC = utilities.c utilities.h

################################################################################

assembler_SOURCES = $(ASSERTS_SRC) $(RUNTIME_SRC) $(UTILITIES_SRC) $(COMMON_SRC) $(SEQUENCES_SRC) $(SLINKLIST_SRC) $(HASHTABLE_SRC) $(CONTIG_SRC) $(GRAPH_SRC) assembler.c

genomewelder_SOURCES = $(ASSERTS_SRC) $(RUNTIME_SRC) $(UTILITIES_SRC) $(COMMON_SRC) $(SEQUENCES_SRC) $(SLINKLIST_SRC) $(HASHTABLE_SRC) $(CONTIG_SRC) genomewelder.c

make_template_SOURCES = $(ASSERTS_SRC) $(RUNTIME_SRC) $(UTILITIES_SRC) $(COMMON_SRC) $(SEQUENCES_SRC) make_template.c

################################################################################

bin_SCRIPTS = fasta.py genomewalker trim_assembly best_hit

################################################################################
