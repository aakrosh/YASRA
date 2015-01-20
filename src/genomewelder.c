/* This module takes as input a fasta file of sequences. The output is another*
 * fasta file of sequences, with the overlapping contigs merged.              *
 * Algorithm:
 * ----------
 * Input : FASTA file with consensus sequence
 * Output: FASTA file with modified consensus sequence
 * Assume: Contigs are ordered as per their hit on the reference sequence.
 *
 * For every pair of adjacent contigs A,B:
 *  align the last "x" bp of A with first "x" bp of B 
 *  if the score of the alignment exceeds a threshold T1, and the pid of the
 *  alignment exceeds a threshold T2, then merge contigs A and B.
 */

# include "utilities.h"
# include "sequences.h"
# include "contig.h"

/* do we want to print the debug information */
int debug_flag = 0;

static float* V = NULL;
static float* F = NULL;
static int* I   = NULL;

static int* II  = NULL;
static int imax = 0;
static int jmax = 0;
static int alignmentscore = 0;

static uint readends = 250;

/*  Align the last 250 bps of the first read and the first 250 bp of the second
 *  read. If the reads are shorter than that, then we can just align the two
 *  reads. 
 */
static void align_read_ends(const seqread* const r1,
                            const seqread* const r2,
                            int* const numaligned, 
                            int* const pid, 
                            int* const continuity)
{
    forceassert(r1 != r2);
    forceassert(r1 != NULL);
    forceassert(r2 != NULL);

    float mscore = 0;

   /* details of the first node */    
    uint cnt1    = slcount(r1->clear);
    uint  l1     = cnt1 > readends ? readends + 1 : cnt1 + 1;
    element* el1 = cnt1 > readends 
                 ? slelement(r1->clear, cnt1 - readends + 1) : r1->clear;

    /* details of the second node */    
    uint cnt2    = slcount(r2->clear);
    uint  l2     = cnt2 > readends ? readends + 1 : cnt2 + 1;
    element* el2;

    V = ckrealloc(V, l1 * l2 * sizeof(float));
    memset(V, 0, l2 * sizeof(float));
    
    I = ckrealloc(I, l1 * l2 * sizeof(int));
    II= ckrealloc(II, l1 * l2 * sizeof(int));
 
    F = ckrealloc(F, l2 * sizeof(float)); 
    memset(F, 0, l2 * sizeof(float));

    float E = 0, G;

    uint i, j;
    uint max1 = 0, max2 = 0;
    for(i = 1; i < l1; i++){

        E = 0;
        V[i * l2] = 0;
        el2 = r2->clear;
        for(j = 1; j < l2; j++){
            /* along the diagonal? */ 
            G = toupper(el1->base) == toupper(el2->base) ? MATCH : MISMATCH; 
            G += V[((i - 1)* l2) + j - 1];
    
            /* gap in read r1? */
            E = (V[(i * l2) + j - 1] - GAPOPEN) > E 
              ?  V[(i * l2) + j - 1] - GAPOPEN - GAPEXTEND : E - GAPEXTEND;

            /* gap in read r2? */
            F[j] = (V[(i - 1) * l2 + j] - GAPOPEN) > F[j] 
                 ?  V[(i - 1) * l2 + j] - GAPOPEN - GAPEXTEND 
                 :  F[j] - GAPEXTEND;
       
            /* which is the path of least resistance? */
            if(G >= MAX(E,F[j])){
                V[(i * l2) + j]  = G;
                I[(i * l2) + j]  = 0;
                II[(i * l2) + j] = 0;
            }else if(E > F[j]){
                assert(E > G);
                V[(i * l2) + j]  = E;
                I[(i * l2) + j]  = 1;
                II[(i * l2) + j] = 1;
            }else{
                V[(i * l2) + j]  = F[j];
                I[(i * l2) + j]  = 2;
                II[(i * l2) + j] = 2;
            }
 
             if(i == l1 - 1 && V[(i * l2) + j] >= mscore){
                mscore = V[(i * l2) + j];   
                max1 = i;
                max2 = j;
            }
            
            el2 = el2->next;
        }
        el1 = el1->next;
    }

    float score = mscore;
    int dir     = 0;
    int flag    = -1;

    i = max1;
    j = max2;
    
    // lets save these values, so we can go through this alignment later.
    imax = max1;
    jmax = max2;
    alignmentscore = mscore;

    int num = 0;
    int mms = 0;
    int gaps= 0;

    while(score > 0){
        dir = I[(i * l2) + j];
        I[(i * l2) + j] = flag;

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            num++;
            if(score < V[((i - 1) * l2) + (j - 1)]) mms++;
            i--; j--;
        }else if(1 == dir){
            gaps++;
            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            gaps++;
            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        flag  = dir;
        score = V[(i * l2) + j];
    }

    *numaligned = num;
    *pid = 0;
    *continuity = 0;
    if(*numaligned > 0){ 
        *pid = (*numaligned - mms) * 100.0 / *numaligned;  
        *continuity = *numaligned * 100.0 / (num + gaps);
    }
}

static int find_fread_end(const seqread* const r1 UNUSED, 
                          const seqread* const r2)
{
    // lets get the alignment between the two reads on here first
    int i = imax;
    int j = jmax;
    float score = alignmentscore;
    int dir     = 0;
    int flag    = -1;

    uint cnt2 = slcount(r2->clear);
    uint ll2  = cnt2 > readends ? readends + 1 : cnt2 + 1;   

    while(score > 0){
        dir = II[(i * ll2) + j];
        II[(i * ll2) + j] = flag;

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            i--; j--;
        }else if(1 == dir){
            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        flag  = dir;
        score = V[(i * ll2) + j];
    }
    return i - 1;
}

// we still have the details of the alignment in V, and the end positions of the
// alignments in imax and jmax. So I create a new read which is essentially a
// composite from the two reads
static seqread* merge_reads(const seqread* const r1, const seqread* const r2)
{
    uint  l1 = slcount(r1->clear) + 1;
    uint  l2 = slcount(r2->clear) + 1;

    element** Elements1 = ckalloc(l1 * sizeof(element*));
    element** Elements2 = ckalloc(l2 * sizeof(element*));

    element* e;
    uint indx;
    for(indx = 1, e = r1->clear; e; indx++, e = e->next) Elements1[indx] = e;
    for(indx = 1, e = r2->clear; e; indx++, e = e->next) Elements2[indx] = e;

    // lets create a dummy read
    seqread* r = new_read("N", NULL, 1, r1->name,
                          FALSE, r1->c, r1->s, r2->e, 0, 1, FALSE, FALSE);

    // lets get the alignment between the two reads on here first
    int i = imax;
    int j = jmax;
    float score = alignmentscore;
    int dir     = 0;
    int flag    = -1;

    uint cnt1 = slcount(r1->clear);
    uint pad1 = cnt1 > readends ? cnt1 - readends : 0;

    element* elements = NULL;

    uint cnt2 = slcount(r2->clear);
    uint ll2  = cnt2 > readends ? readends + 1 : cnt2 + 1;   

    while(score > 0 && i >= 0 && j >= 0){
        dir = II[(i * ll2) + j];
        II[(i * ll2) + j] = flag;
        e = ckalloc(sizeof(element));

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            e->base = Elements1[i+pad1]->base;
            i--; j--;
        }else if(1 == dir){
            e->base = Elements2[j]->base;
            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            e->base = Elements1[i+pad1]->base;
            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        e->qual = 20;
        sladdhead(&elements, e);

        flag  = dir;
        score = V[(i * ll2) + j];
    }
    slreverse(&elements);

    // now string the bases from read from index i to 0 and on the other end
    // from j to the end of the read
    element* es = NULL;    

    for(indx = 1; indx <= i+pad1; indx++){
        e = ckalloc(sizeof(element));
        e->base = Elements1[indx]->base;
        e->qual = Elements1[indx]->qual;
        sladdhead(&es, e);
    }

    e = sllast(elements);
    e->next = es;
    es = elements;

    for(indx = jmax + 1; indx < l2; indx++){
        e = ckalloc(sizeof(element));
        e->base = Elements2[indx]->base;
        e->qual = Elements2[indx]->qual;
        sladdhead(&es, e);
    }

    slreverse(&es);

    ckfree(r->clear);
    r->clear = es;
    
    ckfree(Elements1);
    ckfree(Elements2);
    return r;
}

// free the resources associated with this read and assign the memory as NULL
static void free_read(seqread**  pread)
{
    seqread* r = *pread;
    ckfree(r->name);
    slfreelist(&r->lunused);
    slfreelist(&r->clear);
    slfreelist(&r->runused);
    ckfree(r);
    *pread = NULL;
}

static void print_read(const seqread* const r)
{
    int indx   = 0;
    element* e = r->clear;

    printf(">%s", r->name);
    for(; e; e = e->next, indx++){
        if(indx % 60 == 0){
            printf("\n");
        }
        printf("%c", e->base);
    }
    printf("\n");
}

void welder(const char* const consensusname, 
            const bool iscircular UNUSED,
            const int alignthreshold UNUSED,  
            const int pidthreshold UNUSED, 
            const int conthreshold UNUSED)
{
    // read and save the sequences from this consensusname file
    int numreads = 1;
    seqread** reads = ckalloc(numreads * sizeof(seqread*));
    sequence* sp = read_fasta_sequence(consensusname);

    while(sp){
        reads[numreads - 1] = new_read((char*)sp->sequence, NULL, 
                                       sp->slen, (char*)sp->header, 
                                       FALSE, numreads, -1, -1, 0, sp->slen, 
                                       FALSE, FALSE);      

        sp = get_next_sequence(sp);

        reads = ckrealloc(reads, ++numreads * sizeof(seqread*));
    }
   
    numreads--;
    timestamp("Read all contigs. ");    

    // align the adjacent reads and make an assembly
    int indx  = 0;
    int numaligned = 0;
    int pid = 0; // percent identity in the aligned region, expressed as %
    int continuity = 0; // continuity expressed as a %

    seqread* r = NULL;
    while(indx < (numreads -1)){
        // align the last "x" bases from the first read and the first "x" bases
        // of the second read. If the number of aligned bases is greater than a
        // threshold A, the pid in the aligned region is greater than a
        // threshold P, then find the consensus sequence of the overlap and
        // merge the two sequences
        align_read_ends(reads[indx], reads[indx + 1], 
                        &numaligned, &pid, &continuity);
        if(1 == debug_flag){ 
            fprintf(stderr, "%s %s -> %d %d %d (%d %d)\n", 
            reads[indx]->name,
            reads[indx+1]->name,numaligned,pid,continuity,imax, jmax);
        }

        if(numaligned > alignthreshold && 
           pid        > pidthreshold   && 
           continuity > conthreshold){
            // create a read with the information about the original read and
            // the alignment of their ends.
            r = merge_reads(reads[indx], reads[indx + 1]);
            if(1 == debug_flag){
                fprintf(stderr, "Merge %s and %s\n", 
                reads[indx]->name, reads[indx+1]->name);
            }

            free_read(reads + indx);
            free_read(reads + indx + 1);

            reads[indx + 1] = r;
        }    
       
        indx += 1;
    }

    if(TRUE == iscircular){
        // find the first contig which exists
        int i = 0;
        for(i = 0; i < numreads; i++){
            if(reads[i] != NULL) break;
        }
        if(i == indx) goto done;

        align_read_ends(reads[indx], reads[i], 
                        &numaligned, &pid, &continuity);
        if(1 == debug_flag){
            fprintf(stderr, "%s %s -> %d %d %d\n",
            reads[indx]->name, reads[i]->name, numaligned, pid, continuity);
        }

        if(numaligned > alignthreshold &&
           pid        > pidthreshold   &&
           continuity > conthreshold){
            int j = find_fread_end(reads[indx], reads[i]);

            element* e = reads[indx]->clear;
            int k = 0;
            for(; e && k < j; e = e->next, k++);
            if(1 == debug_flag){
                fprintf(stderr, "Cutting the last contig short\n");
            }
            if(e != NULL) e->next = NULL;
        }

    }

done:
    // print out the details of the updated assembly
    for(indx = 0; indx < numreads; indx++){
        if(reads[indx] == NULL) continue;

        print_read(reads[indx]);
    }

    for(indx = 0; indx < numreads; indx++){
        if(reads[indx] == NULL) continue;
        free_read(reads + indx);
    }
    ckfree(reads);
}

int main(int argc, char** argv)
{
    argv0 = "genomewelder";
    int c;

    bool circulardataset = FALSE;

    int alignthreshold = 10;
    int pidthreshold   = 95;
    int conthreshold   = 90;

    while (1){
        static struct option long_options[] = {
            {"debug"   ,   no_argument      , 0, 'd'},
            {"circular",   no_argument      , 0, 'o'},
            {"numalign",   required_argument, 0, 'n'},
            {"identity",   required_argument, 0, 'i'},
            {"continuity", required_argument, 0, 'c'},
            {0, 0, 0, 0}
        };
  
        int option_index = 0;
        c = getopt_long (argc,argv,"don:i:c:",long_options,&option_index);

        if (c == -1) break;

        switch (c){
            case 0:
                break;
            case 'd':
                debug_flag = 1;
                break;
            case 'n':
                alignthreshold = atoi(optarg);
                break;
            case 'i':
                pidthreshold = atoi(optarg);
                break;
            case 'c':
                conthreshold = atoi(optarg);
                break;
            case 'o':
                circulardataset = TRUE;
                break;
            case '?':
                break;
            default:
                abort();
        }
    }

    /* allocate resources for heavy recursion in this module */
    allocate_resources();

    /* start clock book-keeping */
    t0 = time(0);

    /* the input file with the consensus sequence */
    char* consensusname = argv[optind]; 

    welder(consensusname, circulardataset, 
           alignthreshold, pidthreshold, conthreshold);

    /* this is primarily for Valgrind. Relinquish all the static variables */
    ckfree(V);
    ckfree(F);
    ckfree(I);
    ckfree(II);

   /* print the relevant stats used by the program */
    print_usage();

    return EXIT_SUCCESS;
}
