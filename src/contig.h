#ifndef CONTIG_H
#define CONTIG_H

#include <math.h>

#include "utilities.h"
#include "sequences.h"
#include "slinklist.h"
#include "hashtable.h"

/* the match, mismatch, gapopen and gapextend scores */
#define MATCH 20
#define MISMATCH -20
#define GAPOPEN 80
#define GAPEXTEND 20

/* the threshold score for alignment */
#define THRESHOLD (10 * MATCH)

/* the weights of the two scoring functions used for realignment */
/* 1 only if base is not equal to consensus */
#define WEIGHTSCORE1 0.5
/* score for aligning a symbol x of S with a column is equal to the fraction of
 * symbols not equal to x*/
#define WEIGHTSCORE2 0.5

typedef struct element_st
{
    struct element_st* next;
    struct element_st* prev;
    uchar base;
    uchar qual;
}element;

typedef struct seqread_st
{
    /* name for this read */
    char* name;

    /* where does this align on the reference */
    uint c;       // chromosome 
    uint s;       // start of the alignment
    uint e;       // end of the alignment
    bool complemented;
    
    /* the actual elements of this read */
    element* lunused;
    element* clear;
    element* runused;

    /* where does this read lie on the assembled contig */
    struct contig_st* contig;
    struct column_st* start;
    struct column_st* end;
}seqread;

typedef struct contig_st
{
    struct contig_st* next;
    struct column_st* columns;

    uint index;
    
    uint numreads;
    seqread** reads;  

    /* where does this contig lie on the reference */
    uint c;
    uint s;
    uint e;

    /* this would be set if this contig was marked as "BAD" */
    bool badcontig;
}contig;

typedef struct column_st
{
    struct column_st* next;
    struct column_st* prev;
    
    uint numelements;
    element** elements;
}column;

typedef struct assembly_st
{
    contig* contigs;
}assembly;

/* create a new read */
seqread* new_read(char* sequence,
                  char* quality,
                  const int slen,
                  const char* const name,
                  const bool rc, 
                  const uint c1,
                  const uint s1, const uint e1,
                  const uint s2, const uint e2,
                  const bool doorient,
                  const bool doconvert);

/* create a new contig from this read */
contig* new_contig(assembly* a, seqread* const r);

/* align the two reads and return the score */
float get_read_alignment_score(const seqread* const r1, 
                               const seqread* const r2,
                               int* const pmax1,
                               int* const pmax2,
                               int* const pdiffs);
/* create a new assembly */
assembly* new_assembly();

/* align the reads with each other in context of this assembly */
void align_nodes(assembly* assmbl, 
                 seqread* const r1, 
                 column* r1start,
                 seqread* const r2,
                 const int weight);

/* realign the assembly */
void realign_assembly(assembly* const assmbl,
                      const bool realignonce);

/* mark the contigs that are "BAD" i.e. the ones I should not bother printing
 * out to the user */
void mark_bad_contigs(assembly* const assmbl,
                      const int minavgcov);

/* print the details of this assembly */
void print_assembly(const assembly* const a, 
                    FILE* const fp, 
                    const bool includegaps,
                    hashtable* const map);

/* print the assembly in the ACE format*/
void print_assembly_ace(const assembly* const assmbl, 
                        FILE* const acefile);

/* print the assembly in the SAM format */
void print_assembly_sam(const assembly* const assmbl,
                        FILE* const samfile,
                        hashtable* const map);

/* free the resources from this assembly */
void free_assembly(assembly** pa);

/* free all the static variables in this file. Primarily I dont want to see any
 * more "still reachable" records in valgrind */
void free_alignment_resources();

#endif
