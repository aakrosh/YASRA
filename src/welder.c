/* This module checks to see if any of the gaps in the assembly can be closed.
 * It iterates through pairs of consecutive contigs given by the assembler and
 * merges them if it they satisfy certain creiterion. The -solexa flag should be
 * set if any of the reads used in the assembly are solexa reads*/

#include "util.h"
#include "seq.h"
#include "contig.h"

#define USE "args = ConsensusR [-solexa] > ConsensusRM"

/*maximum length of the overlap to look for*/
#define MAX_OVERLAP 250

/*default number of rows and columns in dp matrix*/
#define DEFAULT_SIZE 1024

/*different score and pid thresholds*/
#define SCORE_THRESHOLD 200
#define SOLEXA_SCORE_THRESHOLD 80
#define PID_THRESHOLD 90
#define SOLEXA_PID_THRESHOLD 80

/*did we use solexa reads*/
bool Solexa = FALSE;

/*is the genome circular*/
bool Circular = FALSE;

/*set of all the contigs which are inputs*/
READ* Reads = NULL;

/*set of contigs which will be outputs*/
CONTIG* Contigs = NULL;

char buffer[DEFAULT_SIZE];

bool same_bases(const char base, const COLUMN* const column)
{
	bool same = TRUE;
	ELEMENT* element = column->head;
	for(; element; element = element->next){
		if(element->el != base){
			same  = FALSE;
			break;
		}
	}

	return same;
}

double calculate_pid(COLUMN* const scolumn, const READ* const read, int** const V, int** const I, const int len1, const int len2, ELEMENT*** const ppelements, COLUMN** const pstart)
{
	int score = 0;
	int i =0, j = 0;
	int max_i = 0, max_j = 0;

	ELEMENT** elements = *ppelements;
	COLUMN* column = scolumn;
	COLUMN* ncolumn = NULL;
	
	/*lets copy the first MAX_OVERLAP elements of the read for easier access*/
	ELEMENT* element = read->unused;
	for(i = 0; element && (i < MAX_OVERLAP); i++){
		elements[i] = element;
		element = element->next;
	}
	
	/*find the max score*/
	for(i = 0; i < len2; i++){
		if(V[i][len1-1] >= score){
			score = V[i][len1-1];
			max_i = i;
		}
	}
	max_j = len1-1;

	i = max_i;
	j = max_j;
	int dir = -1, flag = -1;
	int same = 0;
	int total = 0;

	while(score>0){
		dir = I[i][j];
		I[i][j] = flag;
		if(0 == dir){
			if(same_bases(elements[i-1]->el, column)){
				same++;
			}
			total++;
			i--;j--;
			ncolumn = column;
			column = column->prev;
		}else if(1 == dir){
			j--;
			ncolumn = column;
			column = column->prev;
		}else{
			i--;
		}
		flag = dir;
		score = V[i][j];
	}
	
	if(column == NULL){
		column = ncolumn;
	}
	*pstart = column;

	return (same*100.0/total);
}

/*construct a rough read by taking a contig and taking the first element in each
 * column. This will be used for a rough alignment*/
READ* construct_read(const CONTIG* const contig)
{
	READ* read = ckallocz(sizeof(READ));
	COLUMN* column = contig->columns;
	ELEMENT* el = NULL;

	sprintf(buffer, ">Contig1 0 %d ", contig->count);
	read->header = copy_string(buffer);


	for(column = contig->columns; column && column->head; column = column->next){
		el = ckallocz(sizeof(ELEMENT));
		el->el = column->head->el;
		el->row = &read->row;
		sladd_head(&read->unused, el);
	}
	slreverse(&read->unused);
	read->row = -1;

	return read;
}

int main(int argc, char** argv)
{
	int score_threshold = SCORE_THRESHOLD;
	double pid_threshold = PID_THRESHOLD;

	while(argc > 2){
		argc--;
		if(0 == strncmp(argv[argc],"-solexa", 7)){
			Solexa = TRUE;
			score_threshold = SOLEXA_SCORE_THRESHOLD;
			pid_threshold = SOLEXA_PID_THRESHOLD;
		}
		if(0 == strncmp(argv[argc],"-circular", 9)){
			Circular = TRUE;
		}
	}

	if(argc != 2){
		fatal(USE);
	}

	/*lets read all the contigs and save them*/
	SEQ* sp = NULL;
	READ* read = NULL;
	READ* next = NULL;
	CONTIG* contig = NULL;
	COLUMN* column = NULL;
	COLUMN* start = NULL;
	int max_length = 0;
	int bt = 0;

	if((sp = seq_get(argv[1])) == NULL){
		fatalf("cannot open %s", argv[1]);
	}
	while(sp){
		if(SEQ_LEN(sp) > max_length){
			max_length = SEQ_LEN(sp);
		}
		/*lets make the header compatible with the %s %d %d %s format.*/
		sprintf(buffer, "%s %d %d %s", SEQ_HEAD(sp), bt, bt+SEQ_LEN(sp)-1, SEQ_HEAD(sp));
		bt += SEQ_LEN(sp);
		ckfree(sp->header);
		sp->header = copy_string(buffer);
		read = new_read(sp, NULL);
		sladd_head(&Reads, read);
		if(!seq_read(sp)){
			break;
		}
	}
	max_length += 2;
	seq_close(sp);
	slreverse(&Reads);

	/*lets initialize the score matrix and some other data structures*/
	initialize_scores((Solexa ? SOLEXA_GAP_OPEN :GAP_OPEN),
	                  (Solexa ? SOLEXA_GAP_EXTEND : GAP_EXTEND), 
					   DEFAULT_SIZE);
	/*dp matrices*/
	int size = DEFAULT_SIZE;
	int** V = alloc2D_int(size, size);
	int*  F = ckallocz(size*sizeof(int));;
	int** I = alloc2D_int(size, size);
	ELEMENT** elements = ckalloc(size*sizeof(ELEMENT*));

	int score = 0;
	double pid = 0;
	int count1 = 0, count2 = 0;
	int i = 0;

	read = Reads;
	count1 = slcount(read->unused);
	contig = construct_contig(count1);
	add_contig(contig, read, 0, count1-1);
	read = read->next;
	
	int expected_overlap = 0;
	int max_i = 0, max_j = 0;

	while(contig && read){
		column = contig->last;
		count2 = slcount(read->unused);
		for(i = 0; i < MIN(MAX_OVERLAP-1, contig->count-1); i++){
			column = column->prev;
		}

		score = -1;
		pid = -1;
		
		/*align the read with the present contig under consideration*/
		if ((score = score_contig_ends(column, i+2, read, MIN(count2, MAX_OVERLAP)+1, &V, &F, &I, FALSE, &max_i, &max_j)) > score_threshold &&
			(pid = calculate_pid(contig->last, read, V, I, i+2, MIN(count2,MAX_OVERLAP)+1, &elements,&start)) > pid_threshold){
		   	next = read->next;
			read->next = NULL;
			
			for(column = start, expected_overlap = 0; column != contig->last; column = column->next, expected_overlap++);

			assert(max_i <= MAX_OVERLAP);
			assert(max_j <= MAX_OVERLAP);
			assert(2*expected_overlap < DEFAULT_SIZE);

			align_special(contig, start, contig->last, read, 0, MIN(2*expected_overlap,slcount(read->unused)-1), &V, &F, &I);
			read = next;
#if DEBUG
			fprintf(stderr,"Score:%d Pid:%3.2f\n", score, pid);
#endif
			continue;
		}

#if DEBUG
		fprintf(stderr,"Score:%d Pid:%3.2f\n", score, pid);
#endif
		
		sladd_head(&Contigs, contig);
		next = read->next;
		read->next = NULL;
		count1 = slcount(read->unused);
		contig = construct_contig(count1);
		add_contig(contig, read, 0 , count1-1);
		read = next;
	}

	/*compare the last contig with the first one if this is a circular genome*/
	if(Circular && (slcount(Contigs) >= 1)){
		read = construct_read(sllast(Contigs));
		
		column = contig->last;
		count2 = slcount(read->unused);
		for(i = 0; i < MIN(MAX_OVERLAP-1, contig->count-1); i++){
			column = column->prev;
		}

		/*align the read with the present contig under consideration*/
		if ((score = score_contig_ends(column, i+2, read, MIN(count2, MAX_OVERLAP)+1, &V, &F, &I, FALSE, &max_i, &max_j)) > score_threshold &&
			(pid = calculate_pid(contig->last, read, V, I, i+2, MIN(count2,MAX_OVERLAP)+1, &elements,&start)) > pid_threshold){
			
			for(column = start, expected_overlap = 0; column; column = column->next, expected_overlap++);
			assert(2*expected_overlap < DEFAULT_SIZE);
			
			align_special(contig, start, contig->last, read, 0, MIN(2*expected_overlap,slcount(read->unused)-1), &V, &F, &I);
			slremove(&Contigs, sllast(Contigs));
		}
	
	}

	/*if we have just one contig but we want to check if the ends overlap (if
	 * its circular), and then chop off some parts of the ends so that they are
	 * not repeated*/
	if(Circular && Contigs == NULL){	
		read = construct_read(contig);

		column = contig->last;
		count2 = slcount(read->unused);
		for(i = 0; i < MIN(MAX_OVERLAP-1, contig->count-1); i++){
			column = column->prev;
		}

		/*align the read with the present contig under consideration*/
		if ((score = score_contig_ends(column, i+2, read, MIN(count2, MAX_OVERLAP)+1, &V, &F, &I, FALSE, &max_i, &max_j)) > score_threshold &&
			(pid = calculate_pid(contig->last, read, V, I, i+2, MIN(count2,MAX_OVERLAP)+1, &elements,&start)) > pid_threshold){
#if DEBUG
			fprintf(stderr,"Had to chop off the end of the contig\n");
#endif
			start->next = NULL;				
		}
	}

	sladd_head(&Contigs, contig);
	slreverse(&Contigs);
	
	for(contig = Contigs; contig; contig = contig->next){
		print_aligner_contig(contig);	
	}

	return EXIT_SUCCESS;
}
