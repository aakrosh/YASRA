#ifndef CONTIG_H
#define CONTIG_H

#include "slist.h"
#include "seq.h"

/*extend the contig by this number of bases in one go*/
#define ALLOC 5000

/*when printing a contig we need so many rows. The number is incremented by the
 * same amount*/
#define NUM_ROWS 5

/*default value of quality of a base*/
#define QUALITY 20

/* scores */
#define MATCH 20
#define MISMATCH -20
#define GAPMISMATCH -10

/*gap scores and read length for normal 454 reads*/
#define GAP_OPEN 80
#define GAP_EXTEND 20
#define MAX_SEQ_LENGTH 300


/*the same for the solexa reads*/
#define SOLEXA_GAP_OPEN 80
#define SOLEXA_GAP_EXTEND 20
#define SOLEXA_MAX_SEQ_LENGTH 50

/*the threshold values for the overlaps*/
#define OVERLAP_THRESHOLD 80
#define SOLEXA_CONTAINED_THRESHOLD 100
#define CONTAINED_THRESHOLD(A) (((A)> 100) ? (300) : (200))

/*if the alignments vary more than this, we penalize it*/
#define FUZZ 10 

/*this would define the gap between consecutive reads in the Realigner format*/
#define BANDWIDTH 17

/*size of the dp matrix (multiple of the max_read_length)*/
#define DP_SIZE 2

typedef struct contig
{
	struct contig* next;		/*the next contig*/
	struct read* reads;			/*the reads in this contig*/
	struct column* columns;		/*the positions on this contig*/
	struct column* last;		/*the last column for this contig, used*/
	int bt;						/*begin position of the hit on the template*/
	int et;						/*end position of the hit on the template*/
	int count;					/*number of bases in the contig*/
	int alloced;				/*number of bases allocated for */
}CONTIG;

typedef struct read
{
	struct read* next;			/*the next read in the contig*/
	char* header;				/*the header in the contig*/
	bool complement;			/*is the read complemented*/
	struct column* scolumn;		/*the start column of the aligned read*/
	struct element* unused;		/*trimmed part of the read on the left side*/
	struct element* trimmed;	/*trimmed part of the read on the right side*/
	struct element* head;		/*the first aligned character in the read*/
	int row;					/*the row in the realigner format*/
	/*row is also used as the weight (overlap score) when we are aligning
	 * the contained reads. The value of row is changed to -1 (which is the
	 * default) after that, because a lot of later routines (print_aligner) rely
	 * on that and have checks in place  to make sure that this is true*/
}READ;

#define READ_HEAD(r) 	((r)->header)

typedef struct column
{
	struct column* next;		/*the next column in the contig*/
	struct column* prev;		/*the prev column in the contig*/
	struct element* head;		/*the first element in the column*/
	int counts[6];				/*the counts of the various bases*/
	int index;					/*the index for each column*/
}COLUMN;

typedef struct element
{
	struct element* next;		/*the next element in the same column*/
	struct element* adjacent;	/*the next element of the same read*/
	char el;					/*the base at this position*/
	char qual;					/*the quality value at this position*/
	int* row;					/*which row does this element belong in*/
}ELEMENT;

typedef struct interval
{
	struct interval* next;		/*the next interval in the list*/
	int bt;						/*beginning of the hit on the template */
	int et;						/*end of the hit on the template */
	struct column* start;		/*the corresponding placement on the contig*/
	struct column* end;			/*the last column of alignment of interval */
	struct contig* contig;		/*which contig do I belong to*/
}INTERVAL;

/*initialise the scoring matrices and some of the constants*/
void initialize_scores(const int open, const int extend, const int length);

/*return a new read structure. The default quality value is 20.*/
READ* new_read(const SEQ* const sp, const int* const quals);

/*return what the contained threshold should be*/
int contained_threshold(const int length);

/*take two sequences, align them and return the maximum score*/
int score_reads(const READ* const A, const READ* const B, const bool suffix, int*** const pV, int** const pF, int*** const pI);

/*take a read and an column from a contig and return the score*/
int score_contig_read(const READ* const seq, const COLUMN* const column, const int len1, int*** const pV, int** const pF, int*** const pI, const bool contained);

/*take a read and an column from a contig and return the scores for aligning the
 * right end of the contig and the beginning of the read*/
int score_contig_ends(const COLUMN* const column, const int len1, const READ* const seq, const int len2, int*** const pV, int** const pF, int*** const pI, const bool contained, int* const max_i, int* const max_j);

/*make a new contig of this size*/
CONTIG* construct_contig(const int size);

/* allocate a new contig */
CONTIG* new_contig();

/*free the contigs*/
void free_contigs(CONTIG** pcontigs);

/*free the resources used by the intervals*/
void free_intervals(INTERVAL** pintervals);

/* just add the whole READ to the contig. This is because the READ is the first 
 * read to the contig */
void add_contig(CONTIG* const contig, READ* const node, const int bt, const int et);

/*align the profile with the given sequence */
void align_read(CONTIG* const contig, READ* seq, INTERVAL** const pintervals,  FILE* const rf, int*** const pV, int** const pF, int*** const pI);

/*align the contained read with the existing contig*/
int align_contained_read(CONTIG* const contig, READ* const seq, const INTERVAL* const interval);

/*align the contig and the read. just a wrapper around the align routine*/
int align_special(CONTIG* const contig, COLUMN* const cstart, COLUMN* const cend, READ* const seq, const int rstart, const int rend, int*** const pV, int** const pF, int*** const pI);

/* print out the contig in a form acceptable to Realigner. This routine also
 * sets the row variable of all the reads. That is useful when and
 * if we print the ace file. (Since we want this to work with the realigner)*/
void print_aligner_contig(CONTIG* const contig);

/*Sort the reads with each contig, based on row and position*/
void sort_reads(CONTIG* const contigs);

/*This sets the indices of the columns once for all the contigs so that we 
 * dont have to use slindex for each read when we sort them*/
void assign_indexes(CONTIG* const contigs);

/*print out the partial ace file*/
void print_ace(FILE* const af, const CONTIG* const contigs);

/*free the resources used by the set of intervals*/
void free_intervals(INTERVAL** pintervals);

/*free the resources held by a read*/
void free_read(READ** pread);

/*free the resources used by a list of contigs*/
void free_contigs(CONTIG** pcontigs);

/***********debugging routines***************/
/*print the read*/
void print_read(const READ* const read);

/*print a column*/
void print_column(const COLUMN* const column);

/*print the score matrix*/
void print_scores(int** const V, const int len1, const int len2);

#endif
