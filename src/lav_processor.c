/* lav_processor -- process the blastz/lastz output to select one position of
 * the read if it is possible
 *
 * syntax: lav_processor [-exact/-sig=number]< input.bz > output.bz
 *
 * 
 * This module looks at all the possible placements of a read. It selects one of
 * the placements if
 * 		the alignment score is significantly higher than the remaining
 * 		placements. "significantly" as defined by the used through the
 * 		sig=number
 * 	If several of the placements are still possible (they are not significantly
 * 	different), then one of the placements is selcted randomly.
 * If -exact is specified then only the reads which have multiple placements
 * with the same alignment score are considered and one of the placements is
 * chosen uniformly.
 *
 * This module does not write the complete a-stanza just the s,b,e lines of the
 * a-stanza. However adding that would be trivial and is not required as of now.
 *
 * */

#include "util.h"
#include "slist.h"

#define USE "lav_processor [-exact/-sig=number]< input.bz > output.bz"

#define SCORE_THRESHOLD 1000

#define ENDCURLY "}"

typedef struct read
{
	char* sstanza;	/*the second line in the sstanza for the lav section*/
	int id;			/*contig id for the read*/
	char* hstanza;	/*the second line in the hstanza for the lav section*/
	struct astanza* astanza;/*the a stanza for the positive orientation*/
}READ;

typedef struct astanza
{
	struct astanza* next;
	bool rc;				/*which orientation does this belong to*/
	int score;				
	int b1;
	int e1;
	int b2;
	int e2;
	struct lstanza* lstanza;
}ASTANZA;

typedef struct lstanza
{
	struct lstanza* next;
	char* ugap;
}LSTANZA;

int score_threshold = SCORE_THRESHOLD;
bool only_exacts = FALSE;

char* reference_sstanza = NULL;	
char* reference_hstanza = NULL;

char buffer[1024];

int num_decisions = 0; /*how many repeat decisions did we make?*/
int unambiguous_placements = 0;

/*sort the reads a stanzas based on this */
int compare(const void* const a, const void* const b)
{
	ASTANZA* c = *(ASTANZA**)a;
	ASTANZA* d = *(ASTANZA**)b;
	
	if(c->score == d->score){
		return (d->e1-d->b1)-(c->e1-c->b1);
	}
	
	return d->score - c->score;
}



/*given the orientation and the read, print all the relevant details*/
void print_details(const READ* const read, const ASTANZA* const stanza, const int begin, const int end, const int id, const int rc)
{
	const ASTANZA* astanza = NULL;
	LSTANZA* lstanza = NULL;

	printf("#:lav\n");
	printf("s {\n");
	printf("%s",reference_sstanza);
	if(0 == rc){
		printf("  %s\" %d %d %d %d\n", buffer, begin, end, rc, id);
	}else{
		printf("  %s-\" %d %d %d %d\n", buffer, begin, end, rc, id);
	}
	printf("}\n");
	
	printf("h {\n");
	printf("%s", reference_hstanza);
	if(0 == rc){
		printf("%s\n", read->hstanza);
	}else{
		read->hstanza[strlen(read->hstanza)-1] = '\0';
		printf("%s (reverse complement)\"\n", read->hstanza);
	}
	printf("}\n");
	
	for(astanza = stanza; astanza; astanza = astanza->next){
		if(astanza->rc == rc){
			printf("a {\n");
			printf("  s %d\n", astanza->score);
			printf("  b %d %d\n", astanza->b1, astanza->b2);
			printf("  e %d %d\n", astanza->e1, astanza->e2);
			for(lstanza = astanza->lstanza; lstanza; lstanza = lstanza->next){
				printf("%s", lstanza->ugap);
			}
			printf("}\n");
		}
	}

	printf("x {\n");
	printf("  n 0\n");
	printf("}\n");
}


/*just print the whole read*/
void print_full_read(READ* const read)
{
	int begin, end, id;
	char* ptr = NULL;

	if(sscanf(read->sstanza," %s %d %d %*d %d ", buffer, &begin, &end, &id) != 4){
		fatalf("wrong sstanza:%s", read->sstanza);
	}
	
	ptr = strrchr(buffer,'"');
	if(*(ptr-1) == '-'){
		*(ptr-1) = '\0';
	}else{
		*ptr = '\0';
	}

	ptr = (strstr(read->hstanza,"(reverse"));
	if(ptr){
		*(ptr-1)='"';
		*ptr = '\0';
	}else{
		read->hstanza[strlen(read->hstanza)-1] = '\0';
	}
	
	/*print details of the normal orientation*/
	if(read->astanza){
		print_details(read, read->astanza, begin, end, id, 0);
		print_details(read, read->astanza, begin, end, id, 1);
	}
}

/*free all the resources used by the read*/
void free_read(READ** pread)
{
	READ* read = *pread;
	LSTANZA* lstanza = NULL;
	ASTANZA* astanza = NULL;

	if(read->astanza){
		for(astanza = read->astanza; astanza; astanza = astanza->next){
			for(lstanza = astanza->lstanza; lstanza; lstanza = lstanza->next){
				ckfree(lstanza->ugap);
			}
			slfreelist(&astanza->lstanza);
		}
		slfreelist(&read->astanza);
	}
	ckfree(read->hstanza);
	ckfree(read->sstanza);
	ckfree(read);
}

/*print the read information with the one placement selected*/
void print_select_read(const READ* const read, const int index)
{	
	/*select the astanza to be printed*/
	LSTANZA* lstanza = NULL;
	ASTANZA* astanza = slelement(read->astanza, index);
	
	int begin, end, id;
	char* ptr = NULL;

	if(sscanf(read->sstanza," %s %d %d %*d %d ", buffer, &begin, &end, &id) != 4){
		fatalf("wrong sstanza:%s", read->sstanza);
	}
	
	ptr = strrchr(buffer,'"');
	if(*(ptr-1) == '-'){
		*(ptr-1) = '\0';
	}else{
		*ptr = '\0';
	}

	ptr = (strstr(read->hstanza,"(reverse"));
	if(ptr){
		*(ptr-1)='"';
		*ptr = '\0';
	}else{
		read->hstanza[strlen(read->hstanza)-1] = '\0';
	}
	
	printf("#:lav\n");
	printf("s {\n");
	printf("%s",reference_sstanza);
	if(astanza->rc == FALSE){
		printf("  %s\" %d %d 0 %d\n", buffer, begin, end, id);
	}else{
		printf("  %s-\" %d %d 1 %d\n", buffer, begin, end, id);
	}
	printf("}\n");
	
	printf("h {\n");
	printf("%s", reference_hstanza);
	if(astanza->rc == FALSE){
		printf("%s\n", read->hstanza);
	}else{
		read->hstanza[strlen(read->hstanza)-1] = '\0';
		printf("%s (reverse complement)\"\n", read->hstanza);
	}
	printf("}\n");
	
	printf("a {\n");
	printf("  s %d\n", astanza->score);
	printf("  b %d %d\n", astanza->b1, astanza->b2);
	printf("  e %d %d\n", astanza->e1, astanza->e2);
	for(lstanza = astanza->lstanza; lstanza; lstanza = lstanza->next){
		printf("%s", lstanza->ugap);
	}
	printf("}\n");

	printf("x {\n");
	printf("  n 0\n");
	printf("}\n");
}

/*process the current read to select the position for the read*/
void process_read(READ** const pread)
{
	READ* read = *pread;

	slsort(&read->astanza, compare);
	
	if(read->astanza->next == NULL){
		unambiguous_placements++;
	}

	int score = read->astanza->score;
	int index = 0;
	ASTANZA* astanza = read->astanza;
	ASTANZA* prev = NULL;

	while(astanza && 
	      (score - astanza->score) <= score_threshold){
		index++;
		prev = astanza;
		astanza = astanza->next;
	}

	if(index != 1 &&
	   read->astanza->score != prev->score){
	   num_decisions++;
	}

	if(only_exacts && index != slcount(read->astanza)){
		print_full_read(read);
		return;
	}
	
	/*choose one of the candidates randomly*/
	int select = index * ((double)rand()/(RAND_MAX+1.0));

	print_select_read(read, select);
}

int main(int argc, char** argv)
{
	argv0 = "lav_processor";

	if(argc > 2){
		fatal(USE);
	}
	
	if(argc == 2){
		argc--;
		if(strncmp(argv[argc],"-exact", 6) == 0){
			score_threshold = 0;
			only_exacts = TRUE;
		}
		if(strncmp(argv[argc],"-sig=", 5) == 0){
			score_threshold = atoi(argv[argc]+5);	
		}
	}

	/*getline stuff*/
	unsigned long n = 1;
	char* lineptr = ckalloc(n+1);
	ssize_t num_read = 0;

	/*print verbatim the beginning of the file and the d-stanza*/
	while((num_read = getline(&lineptr, &n, stdin)) != -1){
		printf("%s", lineptr);
		if(strncmp(lineptr, ENDCURLY, 1) == 0){
			break;
		}
	}
	
	/*initialise and seed the random number generator*/
	long seed = time(NULL) + getpid();
	srand(seed);

	/*ignore the first #:lav*/
	if((num_read = getline(&lineptr, &n, stdin)) == -1){
		fatal("incorrect lav file");	
	}

	READ* read  = NULL;	/*which read am i considering*/
	int id = 0;		/*unique identifier of the sequence*/
	int rc = 0;		/*1 = reverse complement*/

	ASTANZA* astanza = NULL;
	LSTANZA* lstanza = NULL;
	
	while(1){
		if((num_read = getline(&lineptr, &n, stdin)) != -1 && 
		   (strncmp(lineptr, "m {", 3) == 0 ||
		    strncmp(lineptr, "#:eof", 5) == 0)){
			break;
		}

		if(strncmp(lineptr,"s {", 3) == 0){
			if(getline(&lineptr, &n, stdin) != -1 && 
			   reference_sstanza == NULL){
			   reference_sstanza = copy_string(lineptr);
			}
			if(getline(&lineptr, &n, stdin) != -1 &&
			   (read == NULL ||
			   (sscanf(lineptr," %*s %*d %*d %d %d ", &rc, &id) == 2 &&
			    id != read->id))){
				/*free the older read*/
				if(read){
					process_read(&read);
					free_read(&read);
				}
				/*allocate a new one*/
			   	read = ckallocz(sizeof(READ));
				read->sstanza = copy_string(lineptr);
				if(sscanf(read->sstanza, " %*s %*d %*d %*d %d", &read->id) != 1){
					fatal("error reading the id for the read");
				}
			}
		}

		if(strncmp(lineptr,"h {", 3) == 0){
			if(getline(&lineptr, &n, stdin) != -1 &&
			   reference_hstanza == NULL){
			   	reference_hstanza = copy_string(lineptr);
			}
			if(getline(&lineptr, &n, stdin) != -1 &&
			   read->hstanza == NULL){
			   	read->hstanza = copy_string(lineptr);
			}
		}

		if(strncmp(lineptr,"a {", 3) == 0){
			astanza = ckallocz(sizeof(ASTANZA));
			if(getline(&lineptr, &n, stdin) != -1 &&
			   sscanf(lineptr, " s %d ", &astanza->score) == 1 &&
			   getline(&lineptr, &n, stdin) != -1 &&
			   sscanf(lineptr, " b %d %d ", &astanza->b1, &astanza->b2) == 2 &&
			   getline(&lineptr, &n, stdin) != -1 && 
			   sscanf(lineptr, " e %d %d ", &astanza->e1, &astanza->e2) == 2){
			}else{
				fatalf("error in reading a stanza %s", lineptr);
			}
			while(1){
				if(getline(&lineptr, &n, stdin) != -1 &&
				   strncmp(lineptr, "}", 1) == 0){
				   	break;
				}
				lstanza = ckalloc(sizeof(LSTANZA));
				lstanza->ugap = copy_string(lineptr);
				sladd_head(&astanza->lstanza, lstanza);
			}
			if(0 == rc){
				astanza->rc = FALSE;
			}else{
				astanza->rc = TRUE;
			}
			slreverse(&astanza->lstanza);
			sladd_head(&read->astanza, astanza);
		}
	}
		
	if(read){
		process_read(&read);
		free_read(&read);
	}

	/*print the mstanza and the eof line*/
	printf("m {\n  n 0\n}\n");
	printf("#:eof\n");
	ckfree(lineptr);
	ckfree(reference_sstanza);
	ckfree(reference_hstanza);
	fprintf(stderr,"Number af unambiguous placements:%d\n", unambiguous_placements);
	fprintf(stderr,"Number of decisions made (excluding unambiguous and exact):%d\n", num_decisions);

	return EXIT_SUCCESS;
}
