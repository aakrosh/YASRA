#include "ace.h"

/*call the consensus based on the bases in the column. IF alleles is set to
 * true, then instead of returning the base with the majority vote return 
 * the IUB nomenclature for the bases in the column.*/
char call_consensus(const COLUMN* const column, const bool alleles)
{
	/*traverse through the bases in the column. Call the base with the majority
	 * support. A,C,T,G get precedence over N,Gap */
	const ELEMENT* element = column->head;
	char base = '*';
	int A=0, C=0, G=0, T=0, Gap=0, N=0, max=0;
	while(element){
		switch(element->el){
			case 'A': 
				if(++A >= max){
					max = A;
					base = 'A';
				}
				break;
			case 'C':
				if(++C >= max){
					max = C;
					base = 'C';
				}
				break;
			case 'G':
				if(++G >= max){
					max = G;
					base = 'G';
				}
				break;
			case 'T':
				if(++T >= max){
					max =T;
					base = 'T';
				}
				break;
			case '*':
				if(++Gap > max){
					max = Gap;
					base = '*';
				}
				break;
			case 'N':
				if(++N > max){
					max = N;
					base = 'N';
				}
				break;
			default:
				break;
		}
		element = element->next;
	}

	if(alleles){
		base = get_alleles(base, A == 0 ? FALSE : TRUE ,C == 0 ? FALSE : TRUE, G == 0 ? FALSE : TRUE ,T == 0 ? FALSE : TRUE);
	}

	return base;
}

/*is the contig formed properly, check that each column in the contig has the
 * correct  consensus call*/
static void check_integrity(const CONTIG* const contig, const char* const consensus)
{
	uint i = 0;
	const COLUMN* column = contig->columns;
	for(i = 0; i < strlen(consensus); i++){
		if(consensus[i] != 'N' &&
		   consensus[i] != '*' &&
		   consensus[i] != call_consensus(column, TRUE)){
			print_column(column);
			fatalf("check the contig on position %d: %c %c", i, consensus[i], call_consensus(column,TRUE));
		}
		column = column->next;
	}
}

/*read the next contig and fill in the requisite datatructures*/
CONTIG* read_contig(FILE* af, char** plineptr, size_t* pn, char** pconsensus, char** pname, struct hash* const hash)
{
	CONTIG* contig = NULL;
	READ* read = NULL;
	ELEMENT* element = NULL;
	char* consensus = NULL;

	ELEMENT* element_iter = NULL;
	ELEMENT* element_iter_next = NULL;

	/*contig*/
	char contig_name[MAX_SEQ_NAME];
	int num_bases, num_reads;

	/*read*/
	char read_name[MAX_SEQ_NAME];
	int reads_handled = 0, offset, read_bases;
	char complement;

	/*random access to the columns of the contig*/
	COLUMN** columns = NULL;

	COLUMN* column_iter = NULL;
	READ* read_iter = NULL;
	int start, end, total;
	char* quals = NULL;
	int i =0, j = 0;

	while(getline(plineptr, pn, af) != -1){
		/*a new contig*/
		if(strncmp(*plineptr, "CO", 2) == 0){
			if(sscanf(*plineptr, "CO %s %d %d 1 U\n", contig_name, &num_bases, &num_reads) != 3){
				fatalf("error in reading the CO line:%s", *plineptr);
			}
			*pname = copy_string(contig_name);

			contig = construct_contig(num_bases);
			if(getline(plineptr, pn, af) != -1){
				consensus = copy_string(*plineptr);
				*pconsensus = consensus;
			}
			assert(((int)strlen(*pconsensus)-1) == num_bases);
			*(consensus + num_bases) = '\0';
			reads_handled = 0;

			columns = ckalloc(num_bases * sizeof(COLUMN*));
			for(i = 0, column_iter = contig->columns; column_iter; column_iter = column_iter->next){
				columns[i++] = column_iter;
			}
			assert(i == num_bases);
		}

		/*quality values for the consensus, right now this is meaningless*/
		if(strncmp(*plineptr, "BQ", 2) == 0 &&
		   getline(plineptr, pn, af) != -1){
		   continue;	
		}

		/*the AF stanza, allocate resources for the reads*/
		if(strncmp(*plineptr, "AF", 2) == 0){
			read = ckallocz(sizeof(READ));
			if(sscanf(*plineptr, "AF %s %c %d\n", read_name, &complement, &offset) == 3){
				read->complement = complement == 'U' ? FALSE : TRUE;
				read->header = copy_string(read_name);
				read->row = offset;
			}
			sladd_head(&contig->reads, read);
			if(++reads_handled == num_reads){
				slreverse(&contig->reads);
				read_iter = contig->reads;
				reads_handled = 0;
			}
			continue;
		}

		/*fill in the details of the reads*/
		if(strncmp(*plineptr,"RD", 2) == 0){
			if(sscanf(*plineptr,"RD %s %d 0 0\n", read_name, &read_bases) != 2){
				fatalf("error in reading the RD line:%s", *plineptr);
			}
			assert(strcmp(read_name, read_iter->header) == 0);
				
			quals = hash_must_find_val(hash, read_name);
			
			if(getline(plineptr, pn, af) != -1){
				assert(((int)strlen(*plineptr)-1) == read_bases);
				/*fill in the unused part on the left*/
				for(i = 0, j = 0; (i < read_bases) && ((*plineptr)[i] >= 'a'); i++){
					element = ckallocz(sizeof(ELEMENT));
					element->el = (*plineptr)[i];
					element->qual = quals[j++]-1;
					sladd_head(&read_iter->unused, element);
				}
				slreverse(&read_iter->unused);

				/*lets read the bases in the read*/
				for(; (i < read_bases) && ((*plineptr)[i] < 'Z'); i++){
					element = ckallocz(sizeof(ELEMENT));
					element->el = (*plineptr)[i];
					if(element->el == GAP){
						element->qual = MIN(quals[j]-1, quals[j-1]-1);
					}else{
						element->qual = quals[j++]-1;
					}
					sladd_head(&read_iter->head, element);
				}
				slreverse(&read_iter->head);
				
				/*lets read the unused element towards the other end*/
				for(; i < read_bases; i++){
					element = ckallocz(sizeof(ELEMENT));
					element->el = (*plineptr)[i];
					element->qual = quals[j++]-1;
					sladd_head(&read_iter->trimmed, element);
				}
				assert(quals[j]=='\0');
				slreverse(&read_iter->trimmed);

				assert(slcount(read_iter->unused)+slcount(read_iter->head)+ slcount(read_iter->trimmed) == read_bases);

				/*now lets put these elements in their place in the contig*/
				column_iter = columns[read_iter->row - 1 + slcount(read_iter->unused)];
				element = NULL;
				element_iter = read_iter->head;
				element_iter_next = element_iter->next;
				while(element_iter && element_iter_next){
					if(element){
						element->adjacent = element_iter;
					}
					element = element_iter;
					sladd_head(&column_iter->head, element_iter);
					column_iter = column_iter->next;
					element_iter = element_iter_next;
					element_iter_next = element_iter_next->next;								}
				element->adjacent = element_iter;
				sladd_head(&column_iter->head, element_iter);
				element_iter->adjacent = NULL;

			}
			ckfree(quals);
			continue;
		}
	
		/*the QA line*/
		if(strncmp(*plineptr,"QA", 2) == 0){
			if(sscanf(*plineptr, "QA %*d %d %d %d\n",&total,&start, &end) != 3){
				fatalf("error in reading the QA line:%s", *plineptr);
			}
			assert((start-1) == slcount(read_iter->unused));
			assert((total-end)  == slcount(read_iter->trimmed));
			read_iter = read_iter->next;

			if(++reads_handled == num_reads){
				check_integrity(contig, consensus);
				return contig;
			}
			continue;
		}

	}
		
	ckfree(columns);
	return contig;
}
