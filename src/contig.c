#include "contig.h"

static char buffer[4096];
static int SCM[128][128];
static int base_int[128];
static char int_base[6];
static int gap_open;
static int gap_extend; 
static int max_read_length;
static int max_dp_size;

void initialize_scores(const int open, const int extend, const int length)
{
	int i = 0, j = 0;
	for(i = 0; i < 128; i++){
		for(j = 0; j < 128; j++){
			SCM[i][j] = (i == j) ?  MATCH : MISMATCH;
		}
	}
	for(i = 0; i < 128; i++){
		SCM['_'][i] = SCM[i]['_'] = GAPMISMATCH;
	}
	SCM['_']['_'] = 0;
	
	/*to convert between bases ACTG_N and 012345*/
	base_int['A'] = 0;
	base_int['C'] = 1;
	base_int['G'] = 2;
	base_int['T'] = 3;
	base_int['N'] = 4;
	base_int['_'] = 5;

	int_base[0] = 'A';
	int_base[1] = 'C';
	int_base[2] = 'G';
	int_base[3] = 'T';
	int_base[4] = 'N';
	int_base[5] = '_';

	gap_open = open;
	gap_extend = extend;
	max_read_length = length;
	max_dp_size = DP_SIZE*length;
}

/*free the resources held by one single read*/
void free_read(READ** pread)
{
	READ* iter = *pread;
	ckfree(iter->header);
	slfreelist(&iter->unused);
	if(iter->trimmed != NULL){
		slfreelist(&iter->trimmed);
	}
}

/*return what the contained threshold should be*/
int contained_threshold(const int length)
{
	if(length < SOLEXA_MAX_SEQ_LENGTH){
		return SOLEXA_CONTAINED_THRESHOLD;
	}

	return CONTAINED_THRESHOLD(length);
}

static ELEMENT* new_element(const char base, const int val)
{
	ELEMENT* element = ckallocz(sizeof(ELEMENT));
	element->el = base;
	element->qual = val;
	return element;
}

/*return a new read structure. The default quality value is 20.*/
READ* new_read(const SEQ* const sp, const int* const quals)
{
	READ* read = ckallocz(sizeof(READ));

	char name[MAX_SEQ_NAME];
	int bt = -1, et = -1;
	ELEMENT* element = NULL;
	
	/*the header of the read*/
	if(sscanf(SEQ_HEAD(sp),"%*s %d %d %s ", &bt, &et, name) != 3){
		fatal("error reading the header for the read");
	}
	sprintf(buffer,"%s %d %d", name, bt, et);
	read->header = copy_string(buffer);

	read->complement = FALSE;
	if(strstr(SEQ_HEAD(sp),"reverse") != NULL){
		read->complement = TRUE;
	}

	/*lets add the elements of the read*/
	for(bt = 0; bt < (int)strlen((char*)SEQ_CHARS(sp)); bt++){
		element = new_element(SEQ_AT(sp,bt), quals == NULL ? 
											QUALITY : quals[bt]);
		element->row = &read->row;									
		sladd_head(&read->unused, element);
	}
	slreverse(&read->unused);

	read->head = NULL;
	read->trimmed = NULL;
	read->scolumn = NULL;

	/*row is set by the print_aligner routine*/
	read->row = -1;

	return read;
}

/* score the base against a profile */
static int score_against_profile(const COLUMN* const A, const char B)
{
	int score = 0;
	int i = 0;
	double temp = 0;
	
	int sum = A->counts[0] + A->counts[1] + A->counts[2] + A->counts[3] + A->counts[4] + A->counts[5]; 
	for(i =0; i < 6; i++){
		temp = ((double)(A->counts[i] * SCM[(int)B][(int)int_base[i]])/(double)sum);
		score += (int)temp;
	}

	return score;
}

/* score two sequences returning the "relavant" score and the dp matrices 
 * involved. Please see the algorithm on Page 244 of "Algorithms on Strings,
 * Trees and Sequences" by Dan Gusfield. The implementation can be optimized
 * later.
 * This routine assumes that all of the required matrices have been initialised
 * to 0, and it doesnt check for that. That responsibility and the resposibility
 * to make sure that the lengths are correct is of the caller.
 * NOTE: If the second sequence is "contained" in the first sequence then
 * the score returned is the score at the end of the second sequence. However
 * if it is a "suffix" then return  the score at the end of the first sequence*/
static int score_sequences(const void* const A, const int len1, const READ* const B, const int len2, const bool contained, const bool is_contig, int*** pV, int**  pF, int*** pI, int* const max_i, int* const max_j)
{
	int i = 0, j = 0;
	int max_score = 0, E = 0, G = 0;
	
	if((j = MAX(len1, len2)) > max_dp_size){
		j = MAX(j,2*max_dp_size);
		fprintf(stderr,"Creating a new %dx%d array\n", j,j);

		free2D_int(pV, max_dp_size);
		*pV = alloc2D_int(j,j);

		*pF = ckalloc(j*sizeof(int));

		free2D_int(pI, max_dp_size);
		*pI = alloc2D_int(j,j);

		max_dp_size = j;
	}

	memset(*pF,0,len1*sizeof(int));

	int* F = *pF;
	int** V = *pV;
	int** I = *pI;

	COLUMN* iter = NULL;		/*use this to iter through the columns*/
	ELEMENT* iter1 = NULL;		/*use this to iter through elements for A*/
	ELEMENT* iter2 = NULL;		/*use this to iter through elements for B*/

	if(is_contig){
		iter = (COLUMN*)A;
	}else{
		iter1 = ((READ*)A)->unused;
	}
	iter2 = ((READ*)B)->unused;
	
	for(i = 1; i < len2; i++){
		E=0;
		if(is_contig){
			iter = (COLUMN*)A;
		}else{
			iter1 = ((READ*)A)->unused;
		}
		for(j = 1; j < len1; j++){
			/* was this a substitution */
			G = is_contig ? score_against_profile(iter, iter2->el) : SCM[(int)iter1->el][(int)iter2->el];
			if(is_contig){
				iter = iter->next;
			}else{
				iter1 = iter1->next;
			}
			G += V[i-1][j-1];
			/* gap in the contig */
			E = (V[i][j-1] - gap_open) > E ? V[i][j-1] - gap_open - gap_extend : E - gap_extend;
			/*gap in the consensus sequence?*/
			F[j] = (V[i-1][j] - gap_open) > F[j] ? V[i-1][j] - gap_open - gap_extend : F[j] - gap_extend;

			/*which is the path of least resistance */
			if(G >= MAX(E,F[j])){
				V[i][j] = G;
				I[i][j] = 0;
			}else if(E >= F[j]){
				V[i][j] = E;
				I[i][j] = 1;
			}else{
				V[i][j] = F[j];
				I[i][j] = 2;
			}

			/* is it the max_score */
			if((contained) && (V[i][j] >= max_score)){
				max_score = V[i][j];
				*max_i = i;
				*max_j = j;
			}
		}
		
		/* is it the max_score */
		if(!contained && (V[i][j-1] > max_score)){
			max_score = V[i][j-1];
			*max_i = i;
			*max_j = j-1;	
		}

		iter2 = iter2->next;
	}
	return max_score;
}

/*take two sequences, align them and return the maximum score*/
int score_reads(const READ* const A, const READ* const B, const bool suffix, int*** const pV, int** const pF, int*** const pI)
{
	int max_score = 0;
	
	int max_i = 0, max_j = 0;
	int len1 = slcount(A->unused)+1;
	int len2 = slcount(B->unused)+1;

	max_score = score_sequences(A, len1, B, len2, !suffix, FALSE, pV, pF, pI, &max_i, &max_j);

	/*if the estimated overlap is very different from that observed by blastz,
	 *the we should give it a low score */
	int b1 = -1, b2 = -1;
	int e1 = -1, e2 = -1;
	sscanf(READ_HEAD(A), "%*s %d %d ", &b1, &e1); 
	sscanf(READ_HEAD(B), "%*s %d %d ", &b2, &e2); 
	if(suffix){
		if((abs((e1-b2) - max_i) > FUZZ)){
			max_score =  0;
		}
	}

	return max_score;
}

/*take a read and an column from a contig and return the score*/
int score_contig_read(const READ* const seq, const COLUMN* const column, const int len1, int*** const pV, int** const pF, int*** const pI, const bool contained)
{
	int max_i = -1, max_j = -1;
	int len2 = slcount(seq->unused) + 1;

	return score_sequences(column, len1, seq, len2, contained, TRUE, pV, pF, pI, &max_i, &max_j);
}

/*take a read and an column from a contig and return the scores for aligning the
 * right end of the contig and the beginning of the read*/
int score_contig_ends(const COLUMN* const column, const int len1, const READ* const seq, const int len2, int*** const pV, int** const pF, int*** const pI, const bool contained, int* const max_i, int* const max_j)
{
	return score_sequences(column, len1, seq, len2, contained, TRUE, pV, pF, pI, max_i, max_j);
}

/*construct a new column for the alignment*/
static inline COLUMN* new_column()
{
	COLUMN* column = ckallocz(sizeof(COLUMN));
	column->index = -1;
	return column;
}

/*make a new contig of this size*/
CONTIG* construct_contig(const int size)
{
	COLUMN* column = NULL;
	CONTIG* contig = ckallocz(sizeof(CONTIG));
	contig->alloced = size;
	int i = 0;
	for(i = 0; i < contig->alloced; i++){
		column = new_column();
		if(contig->columns != NULL){
			contig->columns->prev = column;
		}
		sladd_head(&contig->columns, column);
	}
	contig->last = contig->columns;
	return contig;	
}

/* make a new contig */
CONTIG* new_contig()
{
	return construct_contig(ALLOC);
}

static void free_contig(CONTIG** pcontig)
{
	CONTIG* contig = *pcontig;
	COLUMN* column = NULL;
	READ* read =  NULL;
	for(read = contig->reads; read; read = read->next){
		/*just free the stuff that wouldnt be in the columns*/
		slfreelist(&read->unused);
		slfreelist(&read->trimmed);
		read->head = NULL;
		ckfree(read->header);
	}
	slfreelist(&contig->reads);
	for(column = contig->columns; column; column = column->next){
		slfreelist(&column->head);		
	}
	slfreelist(&contig->columns);
}

/*free the contigs*/
void free_contigs(CONTIG** pcontigs)
{
	CONTIG* contig = NULL;

	for(contig = *pcontigs; contig; contig = contig->next){
		free_contig(&contig);	
	}
	slfreelist(pcontigs);
}

/*free the resources used by the intervals*/
void free_intervals(INTERVAL** pintervals)
{
	slfreelist(pintervals);
}

/* just add the whole READ to the contig. This is because the READ is the first 
 * read to the contig */
void add_contig(CONTIG* const contig, READ* const node, const int bt, const int et)
{
	int length = slcount(node->unused);
	ELEMENT* iter = NULL;
	ELEMENT* element = node->unused;
	COLUMN* column = contig->columns;

	while(element){
		iter = element->next;
		element->next = NULL;
		element->adjacent = iter;
		sladd_head(&column->head, element);
		column->counts[base_int[(int)element->el]]++;
		element = iter;
		contig->last = column;
		column = column->next;
	}
	contig->bt = (contig->count == 0) ?  bt : contig->bt;
	contig->et = et;
	contig->count += length;
	contig->reads = node;

	node->unused = NULL;
	node->head = contig->columns->head;
	node->trimmed = NULL;
}

/*add this element to the column*/
static void add_counts(COLUMN** const pcolumn, ELEMENT** const pelement)
{
	ELEMENT* element = *pelement;
	COLUMN* column = *pcolumn;

	assert(element->adjacent == NULL);
	element->adjacent = element->next;
	element->next = NULL;
	sladd_head(&column->head, element);
	column->counts[base_int[(int)element->el]]++;
}

/*insert element *pelement between *pprev and *pnext, uses the next pointer, so
 * be careful.*/
static void insert_element(ELEMENT** const pprev, ELEMENT** const pelement, ELEMENT** const pnext)
{
	ELEMENT* prev = *pprev;
	ELEMENT* element = *pelement;
	ELEMENT* next = *pnext;

	prev->next = element;
	element->next = next;
}

/*insert a new column between the two columns*/
static void insert_new_column(COLUMN** const pprev, COLUMN** const pnext)
{
	COLUMN* prev = *pprev;
	COLUMN* next = *pnext;
	COLUMN* column = new_column();

	ELEMENT* pelement = prev->head;
	ELEMENT* element = NULL;
	void* tmp = NULL;
	int count = 0;

	while(pelement){
		if(pelement->adjacent != NULL){
			element = new_element('_', 20);
			element->row = pelement->row;
			tmp = pelement->adjacent;
			pelement->adjacent = element;
			element->adjacent = (ELEMENT*)tmp;
			count++;
			sladd_head(&column->head, element);
		}
		pelement = pelement->next;
	}
	
	column->counts[base_int['_']] = count;

	prev->next = column;
	column->next = next;
	next->prev = column;
	column->prev = prev;
}

/* the workhorse like score_sequences but this one aligns the reads and helps
 * in modelling the realigner output and fills in more information about the 
 * reads */
static int align(CONTIG* const contig, COLUMN* const cstart, COLUMN* const cend, READ* seq, const int sstart, const int send, INTERVAL** const pintervals, int*** const pV, int** const pF, int*** const pI)
{
	COLUMN* column = cstart;	
	
	int len1 = 0;
	for(; column != cend; len1++, column = column->next);
	column = cstart;
	
	len1 += 2;
	int len2 = send - sstart + 2;

	int i = 0, j = 0, G = 0, et = -1;
	int max_score = -1;
	int max_i = -1, max_j = -1;
	INTERVAL* interval = NULL;
	if(pintervals){
		interval = *pintervals;
	}

	bool contained = FALSE;
	if(sscanf(READ_HEAD(seq),"%*s %*d %d ", &et) == 1 &&
	   et < contig->et){
		contained = TRUE;
	}

	if((max_score = score_sequences(column, len1, seq, len2, contained, TRUE, pV, pF, pI, &max_i, &max_j)) == 0){
		return FAILURE;
	}

	/* let us check some conditions on the max_score. */
	if(contained && (max_score < contained_threshold(len2-2))){
		return FAILURE;
	}

	/*backtrack fixing the consensus sequence*/
	int dir = 0, flag  = -1;
	int score = max_score;
	i = max_i;
	j = max_j;

	/*find the column at the end of the alignment*/
	for(column = cstart; G < (max_j-1); G++, column = column->next);
	assert(column != NULL);
	COLUMN* last = column;

	/*the read elements are stored for easier access*/
	ELEMENT** read = ckalloc((len2-1)*sizeof(ELEMENT*));
	ELEMENT* element = seq->unused;
	int read_index = 0;
	for(read_index = 0; read_index < (len2-1); read_index++){
		read[read_index] = element;
		element = element->next;
	}

	/*if the read is contained we might trim the edges on the side*/
	if(contained){
		read[i-1]->next = NULL;
		if(i <= (read_index - 1)){
			seq->trimmed = read[i];
		}
	}

	int** V = *pV;
	int** I = *pI;
	
	while(score > 0){
		dir = I[i][j];
		I[i][j] = flag;
		if(0 == dir){
			add_counts(&column, &read[i-1]);
			i--;j--;
		}else if(1 == dir){
			element = new_element('_', 20);
			element->row = &seq->row;
			insert_element(&read[i-1], &element, &(read[i-1]->next));
			add_counts(&column, &element);
			j--;
		}else{
			insert_new_column(&column, &column->next);
			column = column->next;
			contig->count++;
			add_counts(&column, &read[i-1]);
			i--;
		}
		column = column->prev;
		flag = dir;
		score = V[i][j];
	}

	/* i and j both point to the last characters used in the alignment*/
	int k = 0;
	for(k = 0; k < i; k++){
		assert(read[k]->adjacent == NULL);
	}
	
	/*which part of the read did not align*/
	if(i > 0){
		seq->head = read[i-1]->next;
		read[i-1]->next = NULL;
	}else{
		seq->head = seq->unused;
		seq->unused = NULL;
	}
	
	if(interval){
		interval->start = column == NULL ? cstart : column->next;
	}
	/*the start column for the read*/
	seq->scolumn = column == NULL ? cstart : column->next;
	if(contained){
		goto clean;
	}

	/* if it is a suffix we have to do more */
	for(i = max_i, column = last->next;i < read_index; i++){
		add_counts(&column, &read[i]);
		contig->last = column;
		column = column->next;
		contig->count++;
	}
	
	contig->et = et;

clean:
	ckfree(read);

	return SUCCESS;
}

/* align the profile with the given sequence. This routine fills in the
 * information about the interval the read is aligned to. */
void align_read(CONTIG* const contig, READ* seq, INTERVAL** const pintervals, FILE* const rf, int*** const pV, int** const pF, int*** const pI)
{
	int i = 0;
	int read_length = slcount(seq->unused);
	int bt = -1, et = -1;
	INTERVAL* interval = ckalloc(sizeof(INTERVAL));
	sscanf(READ_HEAD(seq),"%*s %d %d ", &bt, &et);
	interval->bt = bt;
	interval->et = et;

	/*is it the first sequence in the seq*/
	if(0 == contig->count){
		interval->start = contig->columns;
		add_contig(contig, seq, bt, et);
		interval->end = contig->last;
		interval->contig = contig;
		sladd_head(pintervals, interval);
		seq->scolumn = contig->columns;
		return;
	}

	/*do i need more memory for the contig*/
	if((contig->count + read_length) > contig->alloced){
		COLUMN* column = NULL;
		COLUMN* columns  = NULL;
		int size = MAX(ALLOC, read_length+1);
		for(i = 0; i < size; i++){
			column = new_column();
			if(columns != NULL){
				columns->prev = column;
			}
			sladd_head(&columns, column);	
		}
		/*add it to the end of the list*/
		column = ((COLUMN*)sllast(contig->columns));
		column->next = columns;
		columns->prev = column;
		contig->alloced += ALLOC;
	}
	
	/* we do not want to align the whole contig sequence. So find the start
	 * for the contig sequence */
	INTERVAL* last =  NULL;
	for(last = *pintervals; last; last = last->next){
		if((bt >= last->bt) && (bt <= last->et)){
			break;
		}
	}

	COLUMN* start = last->start;
	assert(slindex(contig->columns,start) < slindex(contig->columns,contig->last));
	assert(slindex(contig->columns, contig->last)+1 ==  contig->count);

	/*now align the READ with the contig*/
	COLUMN* end = NULL;
	if(et >= contig->et){
		end = contig->last;
	}else{
		/*find one interval which has et > read->et*/
		for(last = *pintervals; last; last = last->next){
			if(last->bt < et && last->et > et){
				break;
			}
		}
		end = last == NULL ? contig->last : last->end;
	}

	if(align(contig, start, end, seq, 0, read_length-1, &interval, pV, pF, pI) == SUCCESS){
		sladd_head(&contig->reads, seq);
		interval->end = end;
	}else{
		if(rf){
			fprintf(rf, "%s\n", READ_HEAD(seq));
		}
		slremove(&contig->reads, seq);
		free_read(&seq);
		interval->start = last->start;
		interval->end = last->end;
	}

	interval->contig = contig;
	sladd_head(pintervals, interval);

#if DEBUG
	printf("%s\n", seq->header);
	print_aligner_contig(contig);
#endif

	return;
}

/*align the contig and the read. just a wrapper around the align routine*/
int align_special(CONTIG* const contig, COLUMN* const cstart, COLUMN* const cend, READ* const seq, const int rstart, const int rend, int*** const pV, int** const pF, int*** const pI)
{	
	int i = 0;
	int read_length = slcount(seq->unused);
	COLUMN* column = NULL;

	/*do i need more memory for the contig*/
	if((contig->count + read_length) > contig->alloced){
		COLUMN* columns  = NULL;
		int size = read_length+1;
		for(i = 0; i < size; i++){
			column = new_column();
			if(columns != NULL){
				columns->prev = column;
			}
			sladd_head(&columns, column);	
		}
		/*add it to the end of the list*/
		column = ((COLUMN*)sllast(contig->columns));
		column->next = columns;
		columns->prev = column;
		contig->alloced += size;
	}

	if((read_length-1) == rend){
		return align(contig, cstart, cend, seq, rstart, rend, NULL, pV, pF, pI);
	}

	/*store them*/
	ELEMENT** read = ckalloc(read_length*sizeof(ELEMENT*));
	ELEMENT* element = seq->unused;
	int read_index = 0;
	for(read_index = 0; read_index < read_length; read_index++){
		read[read_index] = element;
		element = element->next;
	}
	
	if(align(contig, cstart, cend, seq, rstart, rend, NULL, pV, pF, pI) == SUCCESS){
		/*we have some more stuff to do*/
		if(read_length > rend){
			for(i = rend+1,column = contig->last->next; i < read_length; i++){
				add_counts(&column, &read[i]);
				contig->last = column;
				column = column->next;
				contig->count++;
			}
		}
	}

	ckfree(read);	
	return SUCCESS;
}

/* print the output in a format for Realigner. This also cleans up the output
 * on stdout to remove the blockers*/
void print_aligner_contig(CONTIG* const contig)
{
	int* row = ckallocz(NUM_ROWS*sizeof(int));	/*which row is free*/
	int size = NUM_ROWS;						/*size of the 'row' struct*/	
	ELEMENT* element = NULL;
	int i = 0, j = 0;
	int count = 0;
	int num_occupied = 0;
	
	/*this part of the code can be removed if the function is going to be called
	 * just once, but lets do this for now, to make sure that when this is
	 * called, all reads have no row assigned*/
#if DEBUG
	READ* read = contig->reads;
	for(; read; read = read->next){
		read->row = -1;
	}
#endif


	COLUMN* column = contig->columns;
	
	while(column && (count < contig->count)){
		element = column->head;
#if DEBUG
		for(i = 0; i < size; i++){
			printf("%d ", row[i]);
		}
		printf("\n");
#endif
		while(element){
			if(*element->row != -1){
				/*this read has a row allocated*/
				assert(row[*element->row] == -1);
				goto move;
			}

			/* we need to find this read a new row. iterate through the 'row'
			 * structure and find an index where row[index] == 0*/
			for(i = 0; i < size; i++){
				if(row[i] == 0){
					break;
				}
			}
			if((i < size) && row[i] == 0){
				/*allocate this row as the one for this read*/
				*element->row = i;
				row[*element->row] = -1;
				goto move;
			}

			/*we need more cells in the 'row' structure*/
			j = size;
			size += NUM_ROWS;
			row = ckrealloc(row,size*sizeof(int));
			for(i = j; i < size; i++){
				row[i] = 0;
			}
			*element->row = j;
			row[j] = -1;

			assert(*element->row != -1);
move:
			/*is this the last element of the read*/
			if(element->adjacent == NULL){
				row[*element->row] = BANDWIDTH;
			}

			element = element->next;
		}

		num_occupied = 0;
		for(i = 0; i < size; i++){
			if(row[i] == -1){
				num_occupied++;
				continue;
			}
			if(row[i] == BANDWIDTH){
				num_occupied++;
				row[i] -= 1;
				continue;
			}
			if(row[i] == 0){
				continue;
			}
			row[i] -= 1;
		}
		assert(num_occupied == slcount(column->head));
		column = column->next;
		count++;
	}

	column = contig->columns;
	count = 0;
	char* sequence = ckalloc((size+1)*sizeof(char));
	while(column && (count < contig->count)){
		num_occupied = 0;
		memset(sequence, ' ', size);
		sequence[size] = '\0';
		element = column->head;
		while(element){
			sequence[*element->row] = element->el;
			element = element->next;
			num_occupied++;
		}
		if(num_occupied == 0){
			break;
		}
		printf("%s\n", sequence);
		column = column->next;
		count++;
	}
	printf("\n");

	ckfree(row);
	ckfree(sequence);
}

static int compare(const void* const a, const void* const b)
{
	READ* c = *(READ**)a;
	READ* d = *(READ**)b;

	if(c->scolumn->index == d->scolumn->index){
		return c->row - d->row;
	}

	return c->scolumn->index - d->scolumn->index;
}

/*Sort the reads with each contig, based on row and position*/
void sort_reads(CONTIG* const contigs)
{
	CONTIG* contig = NULL;
	for(contig = contigs; contig; contig = contig->next){
		slsort(&contig->reads, compare);
	}
}

/*This sets the indices of the columns once for all the contigs so that we 
 * dont have to use slindex for each read when we sort them*/
void assign_indexes(CONTIG* const contigs)
{
	CONTIG* contig = NULL;
	COLUMN* column = NULL;
	int count = 0;

	for(contig = contigs; contig; contig = contig->next){
		count = 0;
		for(column = contig->columns; column; column = column->next){
			assert(column->index == -1);
			column->index = count++;
		}
	}
}

/*print out the partial ace file*/
void print_ace(FILE* const af, const CONTIG* const contigs)
{
	const CONTIG* contig = NULL;
	static int count = 1;
	const READ* read = NULL;
	const ELEMENT* element = NULL;
	int left = 0, mid = 0, right = 0;
	char name[MAX_SEQ_NAME];

	for(contig = contigs; contig; contig = contig->next){
		fprintf(af,"CO %d %d %d\n", count++, contig->count, slcount(contig->reads));
		read = contig->reads;
		while(read){
			left = 0;
			right = 0;
			mid = 0;
			if(sscanf(read->header, "%s %*d %*d ", name) != 1){
				fatalf("Error reading header of a read:%s", read->header);
			}
			/*is there a comma in the name, if yes we should remove it*/
			if(strchr(name,',') != NULL){
				name[strlen(name)-1] = '\0';
			}
			fprintf(af,"AF %s %c %d\n", name, read->complement ? 'C' : 'U', read->row+1);
			fprintf(af,"RD %s L 0 0\n", name);
			for(element = read->unused; element; element = element->next){
				fprintf(af,"%c", tolower(element->el));
				left++;
			}
			for(element = read->head; element; element = element->adjacent){
				fprintf(af,"%c", element->el == '_' ? '*' : element->el);
				mid++;
			}
			for(element = read->trimmed; element; element = element->next){
				fprintf(af,"%c", tolower(element->el));
				 right++;
			}
			fprintf(af,"\n");
			fprintf(af,"QA 1 %d %d %d\n", left+mid+right, left+1 , left+mid);	
			read = read->next;
		}
	}
}

/*************debugging routines****************/

void print_read(const READ* const read)
{
	fprintf(stderr, "%s\n", read->header);
	ELEMENT* element = read->unused;
	while(element){
		fprintf(stderr,"%c", element->el);
		element = element->next;
	}
	fprintf(stderr,"\n");
	element = read->unused;
	while(element){
		fprintf(stderr,"%d ", element->qual);
		element = element->next;
	}
	fprintf(stderr,"\n");
}

void print_column(const COLUMN* const column)
{
	ELEMENT* element = column->head;
	while(element){
		fprintf(stderr, "%c", element->el);
		element = element->next;
	}
	fprintf(stderr,"\n");
}

void print_scores(int** const V, const int len1, const int len2)
{
	int i = 0, j = 0;
	for(i = 0; i < len2; i++){
		for(j = 0; j < len1; j++){
			printf("%d\t", V[i][j]);
		}
		printf("\n");
	}
}
