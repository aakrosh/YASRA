/* This module is responsible for taking the hits file and producing the
 * assembled contigs
 *
 * syntax: assembler hits.fa [-ace=Contigs.ace] [-repeats=repeats.txt] \
 *         [-rejects=rejects.txt] [-ignore-eid] [-max_length=length] -solexa > MAlign
 *
 * The algorithm can briefly be described as follows:
 * 1) The reads are partitioned into non-contained reads and contained reads.
 * 2) The set of non-contained reads is used to create a graph with the edges
 *    between the reads which satisfy certain overlap creiterion.
 * 3) Prims algortihm is run to form the contigs for the assembly, by finding
 *    the maximum spanning tree.   
 * 4) The contained reads are then aligned after being sorted based on their hit
 *    on the template.
 * 5) The contigs are printed out in a format acceptable to the
 *    realigner/consense and if required a partial ACE file is printed. This is
 *    modifed/added to by consense to form the final ACE file.
 * */

#include "util.h"
#include "seq.h"
#include "hash.h"
#include "contig.h"
#include "graph.h"
#include "config.h"

#define USE "arg: hits.fa [-ace=Contigs.ace] [-repeats=repeats.txt [-rejects=rejects.txt] [-ignore-eid] [-max_length=length] [-solexa]] "


bool Solexa = FALSE;  					/*is this solexa data*/
bool Just_extend = FALSE;				/*just use the non_contained reads*/

struct hash* Ignore = NULL; 			/*the reads which should be ignored*/
/*which of the two names should be used to identify unambiguous placements*/
bool Ignore_eid = FALSE;
int Max_Read_Length = MAX_SEQ_LENGTH; 	/*max length of the read*/

/*a temporary store for the contained reads*/
READ* Contained = NULL;

/*the associated intervals with the non-contained reads*/
INTERVAL* Intervals = NULL;

/*the assembled contigs*/
CONTIG* Contigs = NULL;

/*sort the intervals based on the begin hit on the template*/
int interval_compare(const void* const a, const void* const b)
{
	INTERVAL* c = *(INTERVAL**)a;
	INTERVAL* d = *(INTERVAL**)b;
	
	if(c->bt == d->bt){
		return d->et - c->et;
	}
	return c->bt - d->bt;
}

/*use this to sort the contained reads in decreasing order of scores*/
int read_compare(const void* const a, const void* const b)
{
	READ* c = *(READ**)a;
	READ* d = *(READ**)b;

	return d->row - c->row;
}

/*score the list of contained reads and then sort them in decreasing order of
 * scores. */
void score_contained_reads(READ** plist, const INTERVAL* const interval, int*** const pV, int** const pF, int*** const pI)
{
	READ* iter = NULL;
	
	COLUMN* column = NULL;
	COLUMN* start = interval->start;
	COLUMN* end = interval->end;
	int size = 0;
	 
	for(column = start; column != end; size++,column = column->next);
	size += 2;
	
	for(iter = *plist; iter; iter = iter->next){
		iter->row = score_contig_read(iter, interval->start, size, pV, pF, pI, TRUE);
	}

	slsort(plist, read_compare);
	/*lets revert the row value back to -1*/
	for(iter = *plist; iter; iter = iter->next){
		iter->row = -1;
	}
}

void align_contained_reads(CONTIG* const contig, READ** plist, const INTERVAL* const interval, FILE* const rf, int*** const pV, int** const pF, int*** const pI)
{
	READ* iter = *plist;

	while(*plist){
		iter = slpop(plist);
		if(align_special(contig, interval->start, interval->end, iter, 0, slcount(iter->unused)-1, pV, pF, pI) == SUCCESS){
			sladd_head(&interval->contig->reads, iter);
		}else{
			if(rf){
				fprintf(rf,"%s\n", READ_HEAD(iter));
			}
			free_read(&iter);
		}
	}
}

void process_reads(const char* const reads, FILE* const rf, FILE* const af)
{
	SEQ* sp = NULL;
	CONTIG* contig = NULL;
	CONTIG* tcontig = NULL;

	NODE* node = NULL;
	GRAPH* graph = NULL;
	READ* read = NULL;
	int non_contained = 0;
	int num_contained = 0;

	if((sp = seq_get(reads)) == NULL){
		fatalf("cannot open %s", reads);
	}

	if(0 == SEQ_LEN(sp)){
		return;
	}

	/*what is the rejection file*/
	if(rf){
		fprintf(rf,"Reads that were rejected by the assembler:\n");
	}
		
	char name[MAX_SEQ_NAME];
	char igname[MAX_SEQ_NAME];
	int bt, et;
	int max_index = -1;
	int old_et = INT_MAX;
	bool end = FALSE;
	int num = 0; 			/*header for the contig*/
	char* ptr =  NULL;

	READ* holder = NULL;	
	READ* contained = NULL;	
	INTERVAL* interval = NULL;

	/*allocate memory for the dp matrices*/
	int** V = alloc2D_int(DP_SIZE*Max_Read_Length, DP_SIZE*Max_Read_Length);
	int*  F = ckalloc(DP_SIZE*Max_Read_Length*sizeof(int));
	int** I = alloc2D_int(DP_SIZE*Max_Read_Length, DP_SIZE*Max_Read_Length);
	
	while(sp){
		old_et = INT_MAX;
		graph = new_graph();

		while(sp &&
		      sscanf(SEQ_HEAD(sp), "%s %d %d %s ", igname, &bt, &et, name) == 4 &&
			  bt < old_et){	
			
			/*is this read to be ignored*/
			if(Ignore_eid){
				if((ptr = strchr(name,',')) != NULL){
					*ptr='\0';
				}
				strcpy(igname, name);
			}
			if((Ignore != NULL) && 
		   	   (hash_lookup(Ignore, igname) != NULL)){
			}else if(et <= max_index){      	
				/*is it a contained read*/
				num_contained++;
				if(!Just_extend){
					read = new_read(sp, NULL);
					sladd_head(&Contained, read);
				}
			}else{  							
				/*add this to the graph*/
				non_contained++;
				read = new_read(sp, NULL);
				node = add_node(graph, read, read->header);
				process_node(graph, node, TRUE, &V, &F, &I);

				max_index = et;
				old_et = et;
			}

			if(!seq_read(sp)){
				end = TRUE;
				break;
			}
		}

		if(graph->node_list != NULL){
			/*find the contigs for these reads*/
			Contigs = MST(graph, sllast(graph->node_list), &Intervals, rf, &V, &F, &I);
			
			/*let us orient the contigs, intervals and the contained reads*/
			slreverse(&Contigs);
			slreverse(&Contained);
			slsort(&Intervals, interval_compare);
	
			/*now let us align the contained reads to the contigs*/
			if(!Just_extend){
				/*the begin and end for a non-contained read*/
				holder = NULL;
				contained = NULL;
				interval = Intervals;
			

				for(read = Contained, contig = Contigs; contig; contig = contig->next){
					/*is this contained read for this contig*/
					while(read &&
			      		sscanf(READ_HEAD(read),"%*s %d %d ", &bt, &et) == 2 &&
				 		bt >= contig->bt &&
				  		et <= contig->et){
						holder = NULL;  
						contained = NULL;
						/*is this contained read for this interval*/
						while(read &&
			      	  		sscanf(READ_HEAD(read),"%*s %d %d ", &bt, &et) == 2 &&
					  		bt >= interval->bt &&
					  		et <= interval->et){
							/*get this read and put it on a list*/
							contained = slpop(&read); 
							contained->next = NULL;
							sladd_head(&holder, contained);
						}

						if(contained && holder){
							/*score all the reads for this interval*/
							score_contained_reads(&holder, interval, &V, &F, &I);
							
							/*align reads for this interval decreasing order*/
							if(slindex(contig->columns, interval->start) == -1){
								tcontig = contig;
								contig = interval->contig;
								align_contained_reads(contig, &holder,interval, rf, &V, &F, &I);
								contig = tcontig;
							}else{
								align_contained_reads(contig, &holder,interval, rf, &V, &F, &I);
							}

							/*we dont need the holder*/
							slfreelist(&holder);
							assert(holder == NULL);
						}

						interval = interval->next;
					}
				}
			}
			
			for(contig = Contigs; contig; contig = contig->next){
				fprintf(stderr,"Contig%d: %d-%d\n", ++num, contig->bt, contig->et);
				print_aligner_contig(contig);
			}
			
			if(af){	
				assign_indexes(Contigs);
				sort_reads(Contigs);
				print_ace(af, Contigs);
			}
		}
	
		/*free all the resources used*/
		free_graph(&graph);
		free_contigs(&Contigs);
		free_intervals(&Intervals);
		Contained = NULL;

		if(end){
			break;
		}
	};

	/*lets print some stats*/
	fprintf(stderr,"Number of non-contained reads:%d\n", non_contained);
	fprintf(stderr,"Number of contained reads:%d\n", num_contained);

	seq_close(sp);	
	if(Ignore != NULL){
		free_hash(&Ignore);
	}
	ckfree(F);
}

int main(int argc, char** argv)
{
	argv0 = "assembler";
	bool print_ace = FALSE;
	bool ignore  = FALSE;
	bool rejects = FALSE;
	bool given_max = FALSE;
	FILE* af = NULL;
	FILE* rf = NULL;
	FILE* igf = NULL;

	while(argc > 2){
		argc--;
		if(0 == strncmp(argv[argc],"-ace=", 5)){
			print_ace = TRUE;
			af = ckopen(argv[argc]+5, "w");
		}
		if(0 == strncmp(argv[argc],"-repeats=",9)){
			ignore = TRUE;
			igf = ckopen(argv[argc]+9, "r");
		}
		if(0 == strncmp(argv[argc], "-rejects=",9)){
			rejects = TRUE;
			rf = ckopen(argv[argc]+9, "w");
		}
		if(0 == strncmp(argv[argc], "-solexa", 7)){
			Solexa = TRUE;
		}
		if(0 == strncmp(argv[argc], "-extend", 7)){
			Just_extend = TRUE;
		}
		if(0 == strncmp(argv[argc], "-ignore-eid", 11)){
			Ignore_eid = TRUE;
		}
		if(0 == strncmp(argv[argc], "-max_length=", 12)){
			given_max = TRUE;
			Max_Read_Length = atoi(argv[argc]+12);
		}
	}

	if(argc != 2){
		fatal(USE);
	}
	
	/*lets initialize the score and the constants*/
	if(Solexa && !given_max){
		Max_Read_Length = SOLEXA_MAX_SEQ_LENGTH;
	}
	fprintf(stderr,"Expected maximum length of a read: %d\n", Max_Read_Length);
	initialize_scores((Solexa ? SOLEXA_GAP_OPEN :GAP_OPEN),
	                  (Solexa ? SOLEXA_GAP_EXTEND : GAP_EXTEND), Max_Read_Length);


	/*do we have a list of reads we want to throw out?, put them in a hash*/
	if(ignore){
		Ignore = new_hash(8);
		unsigned long n = 1;
		char* lineptr = ckalloc(n+1);
		ssize_t num_read = -1;
		while((num_read = getline(&lineptr, &n, igf)) != -1){
			/*strip the newline*/
			if(lineptr[strlen(lineptr)-1] == '\n'){
				lineptr[strlen(lineptr)-1] = '\0';
			}
			hash_add(Ignore, lineptr, NULL);
		}
		ckfree(lineptr);
	}

	/*go through the reads and make the longest contigs from them*/
	process_reads(argv[1], rf, af);

	if(ignore) fclose(igf);	
	if(print_ace) fclose(af);
	if(rejects) fclose(rf);
	return EXIT_SUCCESS;
}
