// fix mono-nucleotide runs in 454-sequenced sequences.
// Note that both the positions in the ambiguous file and the substitution
//  file are 1 based and not 0

#include "util.h"
#include "seq.h"
#include "hash.h"

#define USE "fix_asm Assembly subsFile AmbiguityFile string_mononuclear [-H -S -E]"

#define HOMOPOLYMER_RUN 4

#define LINE 1024

int base_int[128];

char buf[LINE];

typedef struct read
{
	struct read* next;		/*next in the list*/
	char* name;				/*unique for a contig, part of the header*/
	SEQ* sequence;			/*sequence for the contig*/
}READ;

/* allocate a new read with the sequence given by the parameter */
READ* alloc_read(const SEQ* const seq)
{
	READ* read = ckallocz(sizeof(READ));
	read->sequence = seq_copy(seq);
	sscanf(read->sequence->header, "%s %*d %*d %*s", buf);
	read->name = copy_string(buf);
	return read;
}

/* free a list of reads and set it to NULL*/
void free_read(READ** read)
{
	ckfree((*read)->name);
	seq_close((*read)->sequence);
}


typedef struct ambiguous
{
	char contig[24];		/* contig name */
	int position;			/* position in the contig */
	char base;				/* the base called for that position */
	int counts[5];			/* the counts for various bases at that position*/
} AMB;

struct hash* amb = NULL;

inline char* get_name(char* const header)
{
	char* name = NULL;
	char* p = NULL;

	if((p = strchr(header,'.')) == NULL){
		fatal("Incorrect format .");
	}
	name = copy_string(p+1);
	if((p = strchr(name,':')) == NULL){
		fatal("Incorrect format :");
	}
		*p='\0';
		return name;
}

void read_amb(const char* const amb_file)
{
	AMB* base = NULL;
	FILE* f = ckopen(amb_file, "r");
	char nbuf[LINE];

	/*ignore the first two lines*/
	if(fgets(buf,LINE,f) == NULL ||
	   fgets(buf, LINE, f) == NULL){
	   	fatal("error in amb file");
   }

	while(fgets(buf, LINE, f)){
		base = ckallocz(sizeof(AMB));
		sscanf(buf,"%s %d %c %d %d %d %d %d", base->contig, &base->position, &base->base, &base->counts[0], &base->counts[1], &base->counts[2], &base->counts[3], &base->counts[4]);
		sprintf(nbuf,"%s.%d", base->contig, base->position);
		hash_add(amb, copy_string(nbuf), base);
	};
	fclose(f);
}

inline void encode()
{
	base_int['A']=0;
	base_int['C']=1;
	base_int['G']=2;
	base_int['T']=3;
	base_int['_']=4;
}

AMB* ambiguous_base(const char* const name, const int position)
{
		char nbuf[LINE];
		AMB* result = NULL;
		
		sprintf(nbuf,"%s.%d", name, position);
		if((result = hash_find_val(amb, nbuf)) != NULL){
			hash_remove_el(amb, nbuf);
			return result;
		}
		
		return NULL;
}

char move_next(SEQ* const sp, READ** pnode, char* const x, char* const y, char*  pname, int* const index)
{
	char* contig_name = NULL;	/*name of the contig from the subs file*/
	char* p = NULL; 			/*temp holder for strchr and strstr*/
	int pos = 0;				/*position in assm*/
	char base;					/*base at pos*/
	
	READ* node = *pnode;
	
	if((p = strchr(x,':')) == NULL){
		fatalf("colon %s", x);
	}

	pos = atoi(p+1);
	base = y[0];
		
	contig_name = get_name(x);
	/*move to the correct contig in our assembly*/
	while(!same_string(contig_name, pname)){
		if(0 == *index){
			printf("%s\n", node->sequence->header);
		}
		for(;*index < SEQ_LEN(node->sequence); ++(*index)){
			putchar(node->sequence->seq[*index]);
		}
		putchar('\n');
				
		if(!seq_read(sp)){
			fatalf("Missing required contig %s", pname);
		}
		free_read(pnode);
		*pnode = alloc_read(sp);
		node = *pnode;
		pname = node->name+1;
		*index = 0;
	}
	
	if(0 == *index){
		printf("%s\n", node->sequence->header);
	}
	for(; *index < pos-1; ++(*index)){
		putchar(node->sequence->seq[*index]);
	}

	fflush(stdout);
	return base;
}

char get_complex_base(const bool A, const bool C, const bool G, const bool T)
{	
	if(A){
		if(C){
			if(G){
				if(T){
					return 'N';
				}else{
					return 'V';
				}	
			}else{
				if(T){
					return 'H';
				}else{
					return 'M';
				}
			}
		}else{
			if(G){
				if(T){
					return 'D';
				}else{
					return 'R';
				}
			}else{
				if(T){
					return 'W';
				}else{
					fatal("Simple base A: should not be here");
				}
			}
		}	
	}else{
		if(C){
			if(G){
				if(T){
					return 'B';
				}else{
					return 'S';
				}	
			}else{
				if(T){
					return 'Y';
				}else{
					fatal("Simple base C: should not be here");
				}
			}
		}else{
			if(G){
				if(T){
					return 'K';
				}else{
					fatal("Simple base G: should not be here");
				}
			}else{
				if(T){
					fatal("Simple base T: should not be here");
				}else{
					fatal("Error: should not be here");
				}
			}
		}
	}
	/*should never reach here*/
	return 'N';	
}

int main(int argc, char **argv) 
{
	bool homopolymer = FALSE;
	bool subst = FALSE;
	bool extra_bases =  FALSE;
	
	while(argc > 5){
		--argc;
		if(strncmp(argv[argc], "-H", 2) == 0){
			homopolymer = TRUE;
		}else if(strncmp(argv[argc] ,"-S", 2) == 0){
			subst = TRUE;
		}else if(strncmp(argv[argc], "-E", 2) == 0){
			extra_bases = TRUE;
		}else{
			fatalf("bad command line argument: %s", argv[argc]);
		}
	}

	/* are we using it correctly ? */
	if(argc != 5){
		fatal(USE);
	}

	amb = new_hash(2);

	/*let us encode for the various bases*/
	encode();
	
	SEQ* sp = NULL;					
	READ* node = NULL;				/*the contig from our assembly*/
	int ocount = 0, ncount = 0; 	/*number of bases in a run in both assm*/
	char x[LINE], y[LINE], z[LINE];
	int index = 0;					/*index in the present contig*/
	char base;						/*base on the assm pos*/
	char tbase;						/*base on the template*/
	AMB* amb_base = NULL;			/*the ambiguous base*/
	bool remove_base = FALSE;	/*do I remove the extra base*/
	int i = 0;	
	
	/* for handling the complext bases*/					
	bool A = FALSE, C = FALSE, G = FALSE, T = FALSE; 	
	char cbase;						/*complex base*/

	/*lets read the bases with low confidence values*/
	read_amb(argv[3]);

	/*this is the string used to check for longer monoruns in our assembly*/
	char* chk_string = argv[4];
	
	FILE* fp = ckopen(argv[2], "r");
	
	/*lets read our assembly*/
	if((sp = seq_get(argv[1])) == NULL){
		fatalf("cannot open %s", argv[1]);
	}

	/* is this file empty */
	if(SEQ_LEN(sp) == 0){
		fatal("Empty Assembly");
	}
	node = alloc_read(sp);
		
	while(fgets(buf, LINE, fp)){
		if(strchr(buf,'>') || strstr(buf,"extra") || (strchr(buf,',') && !strstr(buf,"matches") && !strstr(buf,"label"))){
			/*shorter mononuclear run*/
			if(homopolymer && strchr(buf,',')){
				if(sscanf(buf,"%*s %d %*s %s %d %s", &ocount, x, &ncount, y) != 4){
					fatalf("choked on %s", buf);
				}

				if(ocount < HOMOPOLYMER_RUN){
					continue;
				}

				base = move_next(sp, &node, x, y, node->name+1, &index);
				putchar(node->sequence->seq[index++]);

				putchar('\n');
				while(ocount-- > ncount){
					putchar(base);
				}
				putchar('\n');
				continue;
			}

			/*substitution*/
			if(subst && strchr(buf,'>')){
				if(sscanf(buf,"%*s %s %*s %s %s", x,y,z) != 3){
					fatalf("choked on %s", buf);
				}

				base = move_next(sp, &node, z, y, node->name+1, &index);
				if(node->sequence->seq[index] != y[0]){
					fatalf("substitution: difference in assembly and fix at %d: %c vs %c", index, node->sequence->seq[index], y[0]);
				}

				i = index;
				while(node->sequence->seq[i] == base){
					i--;
				}
				i++;
				while(node->sequence->seq[i] == base){
					if((amb_base = ambiguous_base(node->name+1, i+1)) != NULL){
						remove_base = TRUE;
						break;
					}
					i++;
				}
				
				if(remove_base){
					tbase = x[0];
					/*is it a choice between two simple bases or between 
					 *complex one and a simple one*/
					if((tbase == 'A') || (tbase == 'C') || (tbase == 'G') || (tbase == 'T')){				
						if(amb_base->counts[base_int[(int)amb_base->base]] == amb_base->counts[base_int[(int)tbase]]){
							putchar('\n');
							putchar(tbase);
							putchar('\n');
						}else{
							putchar(base);
						}
					}else{
						ocount = amb_base->counts[base_int[(int)amb_base->base]];
						A = (amb_base->counts[0] == ocount);
						C = (amb_base->counts[1] == ocount);
						G = (amb_base->counts[2] == ocount);
						T = (amb_base->counts[3] == ocount);
						cbase = get_complex_base(A,C,G,T);
						putchar('\n');
						putchar(cbase);
						putchar('\n');
					}

					++index;
					remove_base = FALSE;
				}
				continue;
			}

			/*longer mononuclear run*/
			if(extra_bases && strstr(buf, "extra") && strstr(buf,chk_string)){
				if(sscanf(buf, "%s %*s %s", x, y) != 2){
					fatalf("choked on %s", buf);
				}
	
				base = move_next(sp, &node, x, y, node->name+1, &index);
				if(node->sequence->seq[index] != y[0]){
					fatalf("longer: difference in assmebly and fix at %d: %c vs %c", index, node->sequence->seq[index], y[0]);
				}
			
				i = index; ocount = 0;
				while(node->sequence->seq[i] == base){
					i--;
					ocount++;
				}
				i++;
				
				while(node->sequence->seq[i] == base){
					if(ambiguous_base(node->name+1, i+1)){
						remove_base = TRUE;
						break;
					}
					i++;
				}

				if(remove_base){
					putchar('\n');
					++index;
					remove_base = FALSE;
				}
				continue;
			}
		}
	}

	if(0 == index){
		printf("%s\n", SEQ_HEAD(node->sequence));
	}
	
	for(; index < SEQ_LEN(node->sequence); ++index){
		putchar(node->sequence->seq[index]);
	}
	putchar('\n');
	while(seq_read(sp)){
		free_read(&node);
		node = alloc_read(sp);
		index = 0;
		printf("%s\n", SEQ_HEAD(node->sequence));
		for(; index < SEQ_LEN(node->sequence); ++index){
			putchar(node->sequence->seq[index]);
		}
		putchar('\n');
	}
	
	seq_close(sp);
	fclose(fp);

	return EXIT_SUCCESS;
}
