/*consense -- process the realigner/assembler output to find the consensus
 *
 * syntax: consense [-profile=profile.txt] [-threshold=0.6] \
 *    		[-ace=Contigs.ace] [-amb=file.txt] [-alleles]< MAlign > Consensus
 *
 * This module reads in the output of the realigner/assembler and finds the
 * consensus sequence. If the optional flag -profile is passed it also generates
 * a base profile for the contigs, with details of coverage and support for
 * various bases A,C,G,T,_ at each of the base positions and writes it to
 * profile.txt. 
 * If a file is specified with  the option -amb, then the positions where 
 * the consensus base is supported by less than a particular % of the reads
 * covering that position, are listed in that file. The default value of the
 * threshold is 60% but a user-defined value can be passed by using the 
 * -threshold option. 
 *  If -ace option is set the consense also reads the partial ACE file and write
 *  down a complete ACE file which can then be viewed with Hawkeye or Consed.
 *  If -alleles is specified then the consensus call is the IUPAC name of the
 *	degenerate nucleotide in the column.
 * */

#include "util.h"
#include "iupac.h"
#include "charvec.h"
#include "slist.h"

#define USE "consense [-profile=profile.txt] [-threshold=0.6] \
 			[-amb=file.txt] [-ace=Contigs.ace] < MAlign > Consensus"

#define AMB_THRESHOLD 0.6

#define CHUNK 2000

typedef struct contig
{
	struct contig* next;	/*the next contig in the list*/
	char** seqs;			/*the probable candidates at each position*/
	char* seq;				/*the consensus sequence*/
	int num_index;			/*the unique index of the contig*/
	int num_bases;			/*the number of bases in the contig*/
	int num_reads;			/*the number of reads in the contig*/
	int num_alloced;		/*number of bases alloced for*/

	/*use these two to iterate through the contig*/
	int position;			
	int row;
}CONTIG;

typedef struct read
{
	struct read* next;		/*the next read in the list*/
	char* name;				/*unique name of the read*/
	char* left;				/*the trimmed part of the read on the left*/
	char* seq;				/*the aligned part of the read*/
	char* right;			/*the trimmed part of the read on the right*/
	bool complement;		/*is the read complemented*/
	int position;			/*1-based position on the read*/
	int row;				/*1 based row in the seqs*/
	int* index;				/*which contig do I belong to*/
}READ;

char buffer[1024];

/* return the sequence of the next read and update the read->position to put in
 * the correct position of the read*/
char* find_next_read(CONTIG* const contig, READ** const pread)
{
	int i = 0;
	READ* read = *pread;

	int old_position = 0;
	int position = contig->position;
	int row = contig->row;

	/*find the beginning of the next read*/
	while(1){
		if(0 == position){
			if(row < (int)strlen(contig->seqs[0])){
			   	if(contig->seqs[0][row] != ' '){
			   		break;
			   	}else{
					row++;
			   	}
			}else{
				row = 0;
				position++;
			}
		}else{
			for(i = row; i < (int)strlen(contig->seqs[position]); i++){
				if(contig->seqs[position][i] != ' ' && 
				   ((int)strlen(contig->seqs[position-1]) <= i ||
				   contig->seqs[position-1][i] == ' ' )){
				   row = i;
					break;
				}
			}
			if(i == (int)strlen(contig->seqs[position])){
				row = 0;
				position++;
			}else{
				break;
			}
		}
	}
	old_position = position;

	/*copy the sequence*/
#if DEBUG
	fprintf(stderr,"%s %d %d\n", read->name, position+1, read->row);
#endif
	if(row != (read->row-1)){
		fatalf("Something is wrong for %s, %d (%d vs %d)\n", read->name, position, row, read->row-1);
	}
	read->position = position+1;
	int index = 0;
	while(position != contig->num_bases && 
		  row < (int)strlen(contig->seqs[position]) && 
	      contig->seqs[position][row] != ' '){
		buffer[index++] = contig->seqs[position][row] == '_' ? '*' : contig->seqs[position][row];
		position++;
	}
	buffer[index] = '\0';
	
	/*update the position and row markers for the contig*/
	contig->position = old_position;
	contig->row = row+1;

	return copy_string(buffer);
}

/*display the full ace format file*/
void display_ace(FILE* const ace, CONTIG* const contigs, READ* const reads)
{
	CONTIG* iter = contigs;
	READ* read = reads;
	READ* tmp = NULL;

	int total_reads = 0;
	int total_bases = 0;
	int total_contigs = 0;
	int i = 0;

	for(iter = contigs; iter; iter = iter->next){
		total_contigs++;
		total_bases += strlen(iter->seq);
		while(read && (*read->index == iter->num_index)){
			/*find the next read sequence and its correct position*/
			read->seq = find_next_read(iter, &read);
			total_reads++;
#if DEBUG
			fprintf(stderr,"%s\n", read->seq);
#endif
			read = read->next;
		}
	}
	read = reads;

	/*lets write the stuff in the ace file*/
	fprintf(ace,"AS %d %d\n\n", total_contigs, total_reads);
	for(iter = contigs; iter; iter = iter->next){
		/*contig information*/
		fprintf(ace, "CO Contig%d %d %d 1 U\n", iter->num_index, iter->num_bases, iter->num_reads);
		for(i = 0; i < iter->num_bases; i++){
			fprintf(ace,"%c", (iter->seq[i] == '_') ? '*' : iter->seq[i]);	
		}
		fprintf(ace, "\n\n");
		fprintf(ace, "BQ\n");
		/*these quality values are not calculated but they are put because*/
		for(i = 0; i < iter->num_bases; i++){
			if(iter->seq[i] != '_'){
				fprintf(ace,"97 ");	
			}
		}
		fprintf(ace, "\n\n");

		/*reads in this contig information*/
		tmp = read;
		while(read && (*read->index == iter->num_index)){
			fprintf(ace,"AF %s %c %d\n", read->name+1, read->complement ? 'C':'U', read->position-(int)strlen(read->left));
			read = read->next;
		}

		fprintf(ace,"\n");
		while(tmp && (*tmp->index == iter->num_index)){
			fprintf(ace, "RD %s %d 0 0\n", tmp->name+1, (int)(strlen(tmp->seq)+strlen(tmp->left)+strlen(tmp->right)));
			for(i = 0; i < (int)strlen(tmp->left); i++){
				fprintf(ace,"%c", tmp->left[i]);
			}
			for(i = 0; i < (int)strlen(tmp->seq); i++){
				fprintf(ace,"%c", (tmp->seq[i] == '_') ? '*' : tmp->seq[i]);
			}
			for(i = 0; i < (int)strlen(tmp->right); i++){
				fprintf(ace,"%c", tmp->right[i]);
			}
			fprintf(ace, "\n\n");
			fprintf(ace, "QA 1 %d %d %d\n\n", (int)(strlen(tmp->seq)+strlen(tmp->left)+strlen(tmp->right)),(int)strlen(tmp->left)+1, (int)(strlen(tmp->left)+strlen(tmp->seq)));
			tmp = tmp->next;
		}
	}
}



int main(int argc, char** argv)
{
	argv0 = "consense";
	bool print_profile = FALSE;	/*should I print the base profile*/
	bool print_amb = FALSE;		/*should I print the ambiguous base details*/
	bool print_ace = FALSE;
	bool alleles = FALSE;
	FILE* pf = NULL;			/*file with the profile information*/
	FILE* af = NULL;			/*file with the information of ambiguous bases*/
	FILE* ace = NULL;
	char* acefile = NULL;
	double confidence_threshold = AMB_THRESHOLD;

	while(argc > 1){
		--argc;
		if(0 == strncmp(argv[argc], "-profile=", 9)){
			print_profile = TRUE;
			pf = ckopen(argv[argc]+9,"w");
		}
		if(0 == strncmp(argv[argc], "-amb=", 5)){
			print_amb = TRUE;
			af = ckopen(argv[argc]+5, "w");
		}
		if(0 == strncmp(argv[argc], "-ace=", 5)){
			print_ace = TRUE;
			acefile = copy_string(argv[argc]+5);
			ace = ckopen(acefile, "r");
		}
		if(0 == strncmp(argv[argc], "-threshold=", 11)){
			confidence_threshold = atof(argv[argc]+11);
		}
		if(0 == strncmp(argv[argc], "-alleles", 8)){
			alleles = TRUE;
		}
	}

	if(argc != 1){
		fatal(USE);
	}
	
	/*getline stuff*/
	unsigned long n = 1;
	char* lineptr = ckalloc(n+1);
	ssize_t num_read = -1;

	char name[MAX_SEQ_NAME];
	CONTIG* contig = NULL;
	CONTIG* contigs = NULL;
	READ* read = NULL;
	READ* reads = NULL;
	char complement;
	unsigned int read_index = 0, right = 0;
	int old_num_alloced = 0;

	/*lets read the details of the contigs and the reads from the ace file*/
	if(print_ace){
		while((num_read = getline(&lineptr, &n, ace)) != -1){
			/*contig information*/
			if(0 == strncmp(lineptr, "CO", 2)){
				contig = ckallocz(sizeof(CONTIG));
				if(sscanf(lineptr,"CO %d %*d %d ", &contig->num_index, &contig->num_reads) != 2){
					fatalf("error in reading the CO line:%s", lineptr);
				}
				contig->num_alloced  = contig->num_bases;
				contig->seqs = ckallocz(contig->num_alloced*sizeof(char*));
				contig->seq = ckallocz(contig->num_alloced*sizeof(char));

				sladd_head(&contigs, contig);
				continue;
			}
			
			/*read information*/
			if(0 == strncmp(lineptr,"AF", 2)){
				read = ckallocz(sizeof(READ));
				if(sscanf(lineptr, "AF %s %c %d ", name, &complement, &read->row) != 3){
					fatal("error in reading the AF line");
				}
				read->name = copy_string(name);
				if(complement == 'C'){
					read->complement = TRUE;
				}
				read->index = &contigs->num_index;
				sladd_head(&reads, read);
				continue;
			}

			if(0 == strncmp(lineptr, "RD", 2)){
				if(sscanf(lineptr, "RD %s ", name) != 1 ||
				   strcmp(name, read->name) != 0){
				   	fatal("error in reading the RD line");
				}
			 	continue;
			}
		
			if(0 == strncmp(lineptr, "QA", 2)){
				continue;	
			}

			/*scrape the trimmed left part*/
			read_index = 0;
			while(read_index < (strlen(lineptr)-1) && 
			      lineptr[read_index] >= 'a' &&
				  lineptr[read_index] < 'z'){
				buffer[read_index] = lineptr[read_index];
				read_index++;
			}
			buffer[read_index] = '\0';
			read->left = copy_string(buffer);

			while(read_index < (strlen(lineptr)-1) && 
				  lineptr[read_index++] < 'a'){
			};
			right = read_index;

			/*scrape the trimmed right part*/
			while(read_index <= strlen(lineptr) && 
             	  lineptr[read_index-1] >= 'a' &&
				  lineptr[read_index-1] < 'z'){
				buffer[read_index-right]=lineptr[read_index-1];
				read_index++;
			}
			buffer[read_index-right]='\0';
			read->right = copy_string(buffer);
		}
		slreverse(&contigs);
		slreverse(&reads);
		fclose(ace);
	}

	/*lets read the realigned output from stdin*/
	int contig_num = 1, pos = 1, index = 1;
	int A, C, G, T, Gap, N, total = 0;
	int i = 0, max = 0;
	char base = '*';
	char cbase = '*';

	/*print the schema information in the files*/
	if(print_profile){
		fprintf(pf, "Schema:\nPosition in the contig(1 based)\tCalled base\tNumber of reads suppoting 'A'\tNumber of reads supporting 'C'\tNumber of reads supporting 'G'\tNumber of reads supporting 'T'\tNumber of reads supporting '_'\tCoverage\n");
	}
	if(print_amb){
		fprintf(af, "Schema:\nContig\tPosition in the contig(1 based)\tConsensus base\tNumber of reads supporting 'A'\tNumber of reads supporting 'C'\tNumber of reads supporting 'G'\tNumber of reads supporting 'T'\tNumber of reads supporting '_'\n");
	}

	if(print_profile){
		fprintf(pf,"Contig%d\n", contig_num);
	}
	printf(">Contig%d \n", contig_num++);

	/*read till the end of file is reached*/
	if(print_ace){
		contig = contigs;
	}
	while((num_read = getline(&lineptr, &n, stdin)) != -1){
		/*is it the beginnning of the next contig*/
		if((1 == num_read) && (lineptr[0] == '\n')){
			if((num_read = getline(&lineptr, &n, stdin)) != -1){
				if(print_ace){
					contig->num_bases = pos-1;
					contig = contig->next;
				}
				index = 1;
				pos = 1;
				if(print_profile){
					fprintf(pf, "Contig%d \n", contig_num);
				}
				printf("\n>Contig%d \n", contig_num++);
			}else{
				if(print_ace){
					contig->num_bases = pos -1;
				}
				break;
			}
		}

		if(print_ace && (pos > contig->num_alloced)){
			old_num_alloced = contig->num_alloced;
			contig->num_alloced += CHUNK;
			contig->seqs = 	ckrealloc(contig->seqs, contig->num_alloced*sizeof(char*));
			contig->seq = ckrealloc(contig->seq, contig->num_alloced*sizeof(char));
			memset(contig->seqs+old_num_alloced, 0, CHUNK);	
		}

		if(lineptr[strlen(lineptr)-1] == '\n'){
			lineptr[strlen(lineptr)-1] = '\0';
		}
		if(print_ace){
			contig->seqs[pos-1] = copy_string(lineptr);
		}

		A = C = G = T = Gap = N = 0;
		total = 0, max = 0;

		for(i = 0; i < num_read; i++){
			switch(lineptr[i]){
				case 'A' : 
					A++;   
					total++; 
					if(A >= max){
						max = A; 
						base = 'A';
					} 
					break;
				case 'C' : 
					C++;   
					total++; 
					if(C >= max){
						max = C;
						base = 'C';
					}
					break;
				case 'G' : 
					G++;   
					total++; 
					if(G >= max){	
						max = G;
						base = 'G';
					}
					break;
				case 'T' : 
					T++;   
					total++;
					if(T >= max){	
						max = T;
						base = 'T';
					} 
					break;
				case '_' : 
					Gap++; 
					total++; 
					if(Gap > max){		
						max = Gap;
						base = '_';
					}
					break;
				case 'N':
					N++;
					if(N > max){
						max = N;
						base = 'N';
					}
					break;
				default  : 
					break;
			}
		}

		/*consensus base*/
		cbase = alleles ? get_alleles(base, A == 0 ? FALSE : TRUE ,C == 0 ? FALSE : TRUE, G == 0 ? FALSE : TRUE ,T == 0 ? FALSE : TRUE): base;

		if(print_ace){
			contig->seq[pos-1] = cbase;
		}

		/*is the base ambiguous ? */
		if(print_amb && (((double)(max)/(double)(total)) < confidence_threshold)){
			fprintf(af,"Contig%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\n", contig_num-1, index, cbase, A, C, G, T, Gap);
		}
		
		/*whats the concensus?*/
		if(base != '_' ){
			printf("%c", cbase);
			if(print_profile){
				fprintf(pf, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", index, cbase, A, C, G, T, Gap, total+N);
			}
			index++;
		}

		pos++;
	}
	printf("\n");

	if(print_ace){
		ace = ckopen(acefile, "w");
		display_ace(ace, contigs, reads);	
	}

	if(print_profile) fclose(pf);
	if(print_amb) fclose(af); 
	if(print_ace) fclose(ace);			
	return EXIT_SUCCESS;
}
