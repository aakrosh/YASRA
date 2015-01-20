#include "contig.h"

static float* V = NULL;
static float* F = NULL;
static int* I = NULL;
static column** Columns = NULL;
static element** Elements = NULL;

//static void print_read(const seqread* const r)
//{   
//    pre(r != NULL);
//
//    element* e;
//    for(e = r->lunused; e; e = e->next){
//        printf("%c", e->base == '-'? '*' : tolower(e->base));
//    }
//    for(e = r->clear; e; e = e->next){
//        printf("%c", e->base == '-'? '*' : e->base);
//    }
//    for(e = r->runused; e; e = e->next){
//        printf("%c", e->base == '-'? '*' : tolower(e->base));
//    }
//    printf("\n");
//}
//
//static void print_column(const column* const c)
//{
//    uint i;
//    printf("\t%d\t", c->numelements);
//    for(i = 0; i < c->numelements; i++){
//        printf("%c", c->elements[i] == NULL ? 'Z':c->elements[i]->base);
//    }
//    printf("\n");
//}
//
//static void print_contig_columns(const contig* const c)
//{
//    pre(c != NULL);
//
//    column* col;
//    for(col = c->columns; col; col = col->next){
//        assert(col->next == NULL || col->next->prev == col);
//        print_column(col);
//    }
//}

seqread* new_read(char* sequence,
                  char* quality,
                  const int slen,
                  const char* const name,
                  const bool rc, 
                  const uint c1,
                  const uint s1, const uint e1,
                  const uint s2, const uint e2,
                  const bool doorient,
                  const bool doconvert)
{
    pre(sequence != NULL);
    pre(slen > 0);
    pre(name != NULL);

    seqread* r = ckallocz(sizeof(seqread));

    r->name = copy_string(name);

    r->c = c1;
    r->s = s1;
    r->e = e1;
    r->complemented = rc;

    if(doorient == TRUE && rc == TRUE){
        sequence = (char*)reverse_complement_string((uchar*)sequence, slen);    
        if(quality != NULL){
            quality  = reverse_string(quality, slen);
        }
    }

    uint i;
    element* e;

    for(i = 0; i < s2; i++){
        e = ckalloc(sizeof(element));
        e->base = rc == TRUE ? tolower(sequence[i]) : sequence[i];
        if(doconvert == TRUE){
            e->qual = quality == NULL ? 20 : quality[i] - 64;
        }else{
            e->qual = quality == NULL ? 20 : quality[i] - 33;
        }
        sladdhead(&r->lunused, e);
    }
    slreverse(&r->lunused);

    for(i = s2; i < e2; i++){
        e = ckalloc(sizeof(element));
        e->base = rc == TRUE ? tolower(sequence[i]) : sequence[i];
        e->qual = quality == NULL ? 20 : quality[i] - 33;
        sladdhead(&r->clear, e);
    }
    slreverse(&r->clear);

    for(i = e2; i < (uint)slen; i++){
        e = ckalloc(sizeof(element));
        e->base = rc == TRUE ? tolower(sequence[i]) : sequence[i];
        e->qual = quality == NULL ? 20 : quality[i] - 33;
        sladdhead(&r->runused, e);
    }
    slreverse(&r->runused);

    element* p;
    for(p = NULL, e = r->lunused; e; p = e, e = e->next) e->prev = p;
    for(p = NULL, e = r->clear  ; e; p = e, e = e->next) e->prev = p;
    for(p = NULL, e = r->runused; e; p = e, e = e->next) e->prev = p;

    return r;
}

/* align the two reads and return the score */
float get_read_alignment_score(const seqread* const r1, 
                             const seqread* const r2,
                             int* const pmax1,
                             int* const pmax2,
                             int* const pdiffs)
{
    forceassert(r1 != r2);
    forceassert(r1 != NULL);
    forceassert(r2 != NULL);

    float mscore = 0;

   /* details of the first node */    
    uint  l1 = slcount(r1->clear) + 1;
    element* el1;

    /* details of the second node */    
    uint  l2 = slcount(r2->clear) + 1;
    element* el2;

    V = ckrealloc(V, l1 * l2 * sizeof(float));
    memset(V, 0, l2 * sizeof(float));
    
    I = ckrealloc(I, l1 * l2 * sizeof(int));
 
    F = ckrealloc(F, l2 * sizeof(float)); 
    memset(F, 0, l2 * sizeof(float));

    float E = 0, G;

    uint i, j;
    uint max1 = 0, max2 = 0;
    el1 = r1->clear;
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
                V[(i * l2) + j] = G;
                I[(i * l2) + j] = 0;
            }else if(E > F[j]){
                assert(E > G);
                V[(i * l2) + j] = E;
                I[(i * l2) + j] = 1;
            }else{
                V[(i * l2) + j] = F[j];
                I[(i * l2) + j] = 2;
            }
 
             if(V[(i * l2) + j] >= mscore){
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

    int num1 = 0;
    int num2 = 0;
    int diffs= 0;

    while(score > 0){
        dir = I[(i * l2) + j];
        I[(i * l2) + j] = flag;

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            num1++;
            num2++;
            if(score < V[((i - 1) * l2) + (j - 1)]) diffs++;
            i--; j--;
        }else if(1 == dir){
            num1++;
            diffs++;
            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            num2++;
            diffs++;
            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        flag  = dir;
        score = V[(i * l2) + j];
    }

    *pmax1 = num1;
    *pmax2 = num2;
    *pdiffs= diffs;    

    return mscore;
}

/* create a new assembly */
assembly* new_assembly()
{
    assembly* a = ckallocz(sizeof(assembly));
    return a;
}

/* create a new column with the element "e" as the only element */
static column* new_column(element* const e)
{
    pre(e != NULL);

    column* c      = ckallocz(sizeof(column));
    c->numelements = 1;
    c->elements    = ckalloc(sizeof(element*));
    c->elements[0] = e;

    return c;
}

/* create a new contig from this read */
contig* new_contig(assembly* a, seqread* const r)
{
    pre(a != NULL);
    pre(r != NULL);

    static uint index = 1;
    contig* c = ckallocz(sizeof(contig));

    c->index = index++;

    /* the contig covers all the region covered by this read */
    c->c = r->c;
    c->s = r->s;
    c->e = r->e;

    /* add all the columns to the contig */
    column* clns = NULL;
    column* cln;
    element* iter;
    
    for(iter = r->clear; iter; iter = iter->next){
        cln = new_column(iter);
        sladdhead(&clns, cln);
    }
    column* last = clns;
    slreverse(&clns);

    /* add the links to get the prev pointer working */
    column* prev = NULL;
    for(cln = clns; cln; cln = cln->next){
        cln->prev = prev;
        prev = cln;
    }

    c->columns = clns;
    sladdhead(&a->contigs, c);

    /* add the read to the contig */
    c->numreads = 1;
    c->reads    = ckalloc(c->numreads * sizeof(seqread*));
    c->reads[0] = r;

    /* update the read details */    
    r->contig = c;
    r->start  = c->columns;
    r->end    = last;   
    forceassert(r->start != r->end);

    return c;
}

static inline float scorenuc(const column* const el1, const element* const el2)
{
    pre(el1 != NULL);
    pre(el2 != NULL);

    uint i = 0;
    int score = 0;
    for(; i < el1->numelements; i++){   
        if(el1->elements[i] != NULL){
            score = toupper(el1->elements[i]->base) == toupper(el2->base) 
                  ? score + MATCH : score + MISMATCH;
        }
    }
    
    return score * 1.0 / el1->numelements;
}

/* align the reads with each other in context of this assembly */
void align_nodes(assembly* assmbl, 
                 seqread* const r1, 
                 column*  r1start,
                 seqread* const r2,
                 const int weight UNUSED)
{
    forceassert(r1 != NULL);
    forceassert(r2 != NULL);
    forceassert(r1 != r2);
    
    uint i, j, k; /* counters in loops */

    contig* c1 = r1->contig == NULL ? new_contig(assmbl, r1) : r1->contig;
    assert(sldistance(r1->start,r1->end) == slcount(r1->clear));
    contig* c2 UNUSED= NULL;
    forceassert(r2->contig == NULL);  

    /* does this read extend this contig */
    bool extend = FALSE;
    if(r2->e >= c1->e) extend = TRUE;

    /* if the contig has 0 bases from the first read, then it is pretty much
     * pointless to try to align read 2 to this contig */
    if(r1->clear ==  NULL){
        return;
    }

    /*is there more stuff in this contig beyond the last element of this read */
    column* last = r1->end;
    if(extend == TRUE || r2->e >= r1->e){
        last = sllast(last);
    }
    forceassert(last != NULL);

    /* well if this is a contig formed of only one read, then extra bookkeepin*/
    if(r1start == NULL) r1start = r1->start;

    /* details of the read which has a contig and columns for easier access*/
    column* el1;

    uint l1 = sldistance(r1start, last) + 1;
    
    Columns = ckrealloc(Columns, l1 * sizeof(column*));
    for(i = 1, el1 = r1start; el1 && el1 != last; i++, el1 = el1->next){
        Columns[i] = el1;
    }
    forceassert(i == (l1 - 1));
    Columns[i] = el1;
    //change this back to an forceassert
    if(Columns[i] == NULL){ 
        return;
    }
    
    /* details of the read which is to be added */
    uint l2 = slcount(r2->clear) + 1;
    element* el2;

    Elements = ckrealloc(Elements, l2 * sizeof(element*));
    for(i = 1, el2 = r2->clear; el2; i++, el2 = el2->next){
        Elements[i] = el2;
    }
    forceassert(i == l2);

    /* matrices for the dp matrix and variables for the alignments */    
    V = ckrealloc(V, l1 * l2 * sizeof(float));
    memset(V, 0, l2 * sizeof(float));
 
    I = ckrealloc(I, l1 * l2 * sizeof(int));

    F = ckrealloc(F, l2 * sizeof(float)); 
    memset(F, 0, l2 * sizeof(float));

    float E = 0, G;

    /* the maximum score, and the element/column where we get this */
    float mscore = 0;
    uint max1    = 0;
    uint max2    = 0;

    /* the element in read2, beyond which it is not being used on the right */
    element* r2max = NULL;  

    el1 = r1start;
    for(i = 1; i < l1; i++){

        E = 0;
        V[i * l2] = 0;
        el2 = r2->clear;
        for(j = 1; j < l2; j++){
            /* along the diagonal? */ 
            G = scorenuc(el1, el2); 
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
                V[(i * l2) + j] = G;
                I[(i * l2) + j] = 0;
            }else if(E > F[j]){
                forceassert(E > G);
                V[(i * l2) + j] = E;
                I[(i * l2) + j] = 1;
            }else{
                V[(i * l2) + j] = F[j];
                I[(i * l2) + j] = 2;
            }
     
            /* if this read extends the contig then the maximum score has to 
             * be the last element aligned for this read, otherwise find the
             * highest score on a diagonal */
            if(extend == TRUE){
                if(i == (l1 - 1) && V[(i * l2) + j] >= mscore){
                    mscore = V[(i * l2) + j];
                    max1   = i;
                    max2   = j;
                }
            }else{
                if(V[(i * l2) + j] >= mscore){
                    mscore = V[(i * l2) + j];
                    max1   = i;
                    max2   = j;
                    r2max  = el2;
                }
            }
            
            el2 = el2->next;
        }
        el1 = el1->next;
    }

    /* if the two nodes do not align at all, then this read needs its own contig
     */
    if(mscore == 0){
        c2 = new_contig(assmbl, r2);
        return;
    }

    /* align the read with the contig */
    float score = mscore;
    int dir     = 0;
    int flag    = -1;

    i = max1;
    j = max2;

    /* do we need to move some part of the read to the runused section? */
    if(!extend && r2max->next != NULL){
        for(el2 = r2max->next; el2->next; el2 = el2->next);
        el2->next = r2->runused;
        r2->runused = r2max->next;
        r2->runused->prev = NULL;
        r2max->next = NULL;
    }

    while(score > 0){
        dir = I[(i * l2) + j];
        I[(i * l2) + j] = flag;

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            el1 = Columns[i];
            el1->numelements++;
            el1->elements = ckrealloc(el1->elements, 
                                      el1->numelements * sizeof(element*));
            el1->elements[el1->numelements - 1] = Elements[j];
            i--; j--;
        }else if(1 == dir){
            column* prev = Columns[i];

            /* how many of the reads in the previous column do not end there */
            int l = 0;
            for(k = 0; k < prev->numelements; k++){
                if(prev->elements[k]->next != NULL) l++;
            }

            /* a new column for the contig */
            column* col = ckalloc(sizeof(column));
            col->numelements = l + 1;
            col->elements = ckalloc(col->numelements * sizeof(element*));
            for(k = 0, l = 0; k < prev->numelements; k++){
                if(prev->elements[k]->next == NULL) continue;
                element* tmp = prev->elements[k]->next;

                el2 = ckalloc(sizeof(element));
                el2->base = '-';
                el2->qual = MIN(prev->elements[k]->qual,tmp->qual);

                prev->elements[k]->next = el2;
                el2->prev = prev->elements[k];

                el2->next = tmp;
                tmp->prev = el2;
                
                col->elements[l] = el2;
                l++;
            }
            
            el2 = Elements[j];
            col->elements[l] = el2;

            /* add this column in its location */
            column* tmp = prev->next; 
            prev->next = col;
            col->prev = prev;
            col->next = tmp;
            tmp->prev = col;
            
            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            el1 = Columns[i];
            el1->numelements++;
            el1->elements = ckrealloc(el1->elements, 
                                      el1->numelements * sizeof(element*));

            element* tmp = Elements[j]->next;

            /* the new gap element */
            el2 = ckalloc(sizeof(element));
            el2->base = '-';
            el2->qual = Elements[j]->qual;
            if(tmp != NULL){
                el2->qual = MIN(tmp->qual, Elements[j]->qual);
            }
                
            /* insert it in place (in the read) */
            Elements[j]->next = el2;
            el2->prev = Elements[j];
            
            el2->next = tmp;
            if(tmp != NULL) tmp->prev = el2;

            /* add it to the column */
            el1->elements[el1->numelements - 1] = el2;
            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        flag  = dir;
        score = V[(i * l2) + j];
    }

    /* do we need to add to the lunused section of the read? */
    if(j > 0){
        element* tmp = r2->clear;
        for(el2 = r2->clear, k = 0; k < (j - 1); el2 = el2->next, k++);
        r2->clear = el2->next;
        r2->clear->prev = NULL;
        el2->next = NULL;

        if(r2->lunused){
            for(el2 = r2->lunused; el2->next; el2 = el2->next);
            el2->next = tmp;
        }else{
            r2->lunused = tmp;            
        }
    }

    /* update the details for the new read */
    r2->contig = c1;
    forceassert(r2->c == c1->c);
    r2->start  = Columns[i+1];

    /* if this read extends this contig, then we need to add the remaining bases
     * from the read*/
    el1 = NULL;
    j = max2;
    if(extend == TRUE && j < (l2 - 1)){
        while(j < (l2 - 1)){
            column* col = new_column(Elements[++j]);
            sladdhead(&el1, col);
        }
        slreverse(&el1);

        /* settle the back links for this new piece */
        column* prev = NULL;
        column* iter = NULL;
        for(iter = el1; iter; iter = iter->next){
            iter->prev = prev;
            prev = iter;
        }

        forceassert(Columns[max1]->next == NULL);
        Columns[max1]->next = el1;
        el1->prev = Columns[max1];

        /* update the details for the read and the contig */
        r2->end = sllast(el1);
        c1->e = r2->e;
    }else{
        r2->end = Columns[max1];
    }

    forceassert(r1->lunused == NULL || r1->lunused->prev == NULL);
    forceassert(r1->clear->prev == NULL);
    forceassert(r1->runused == NULL || r1->runused->prev == NULL);
    forceassert(r1->start != NULL);
    forceassert(r1->end != NULL);
    forceassert(r2->lunused == NULL || r2->lunused->prev == NULL);
    forceassert(r2->clear->prev == NULL);
    forceassert(r2->runused == NULL || r2->runused->prev == NULL);
    forceassert(r2->start != NULL);
    forceassert(r2->end != NULL);


    /* add this read to this contig */
    c1->numreads++;
    c1->reads = ckrealloc(c1->reads, c1->numreads*sizeof(seqread*));
    c1->reads[c1->numreads - 1] = r2;

    assert(sldistance(r1->start,r1->end) == slcount(r1->clear));
    assert(sldistance(r2->start,r2->end) == slcount(r2->clear));
}

/* Infer the consensus sequence and quality for the contig */
static void infer_consensus(const contig* const contig, 
                            element** const pconsensus)
{
    pre(contig != NULL);

    element* consensus = *pconsensus;
    
    column* col;
    
    /* sum of quality values for the bases in the two strands */
    uint A,C,G,T,N,GAP;
    uint a,c,g,t,n,gap;

    /* counts of bases */
    uint nA, nC, nG, nT, nN, ngap;

    /* the maximum quality value for the base on each strand */
    uint maxA,maxC,maxG,maxT,maxN,maxGAP;
    uint maxa,maxc,maxg,maxt,maxn,maxgap;

    /* base at a location on the read and the quality value at same location */
    char b; 
    uint q; 

    uint i, j;
    element* e;
    for(i = 0, col = contig->columns; col; ++i, col = col->next){
        A = C = G = T = N = GAP = 0;
        a = c = g = t = n = gap = 0;

        nA = nC = nG = nT = nN = ngap = 0;

        maxA = maxC = maxG = maxT = maxN = maxGAP = 0; 
        maxa = maxc = maxg = maxt = maxn = maxgap = 0;

        for(j = 0; j < col->numelements; j++){
            b = col->elements[j]->base;
            q = col->elements[j]->qual;

            switch(b){
                case 'A': 
                    nA += 1;
                    A  += q; 
                    if(q > maxA) maxA = q; 
                    break;
                case 'a':
                    nA += 1; 
                    a  += q;
                    if(q > maxa) maxa = q; 
                    break;
                case 'C':
                    nC += 1; 
                    C  += q;
                    if(q > maxC) maxC = q; 
                    break;
                case 'c':
                    nC += 1;
                    c  += q;
                    if(q > maxc) maxc = q; 
                    break;
                case 'G':
                    nG += 1; 
                    G  += q; 
                    if(q > maxG) maxG = q;
                    break;
                case 'g':
                    nG += 1; 
                    g  += q;
                    if(q > maxg) maxg = q;
                    break;
                case 'T':
                    nT += 1; 
                    T  += q; 
                    if(q > maxT) maxT = q;
                    break;
                case 't':
                    nT += 1; 
                    t  += q;
                    if(q > maxt) maxt = q; 
                    break;
                case 'N':
                    nN += 1; 
                    N  += q;
                    if(q > maxN) maxN = q;
                    break;
                case 'n':
                    nN += 1;
                    n  += q;
                    if(q > maxn) maxn = q; 
                    break;
                case '-':
                    ngap += 1;
                    e = col->elements[j];
                    for(; e && e->base == '-'; e = e->next);
                    if(e == NULL){
                        e = col->elements[j];
                        for(; e && e->base == '-'; e = e->prev);
                    }
                    if(e->base < 'a'){
                        GAP += q;
                        if(q > maxGAP) maxGAP = q;
                    }else{
                        gap += q;
                        if(q > maxgap) maxgap = q; 
                    }
                    break;
                default : fatalf("unknown base: %x\n", b);
            } 
        }
     
        int tA   = (maxA + 0.5 * (A - maxA)) + (maxa + 0.5 * (a - maxa));
        int tC   = (maxC + 0.5 * (C - maxC)) + (maxc + 0.5 * (c - maxc));
        int tG   = (maxG + 0.5 * (G - maxG)) + (maxg + 0.5 * (g - maxg));
        int tT   = (maxT + 0.5 * (T - maxT)) + (maxt + 0.5 * (t - maxt));
        int tN   = (maxN + 0.5 * (N - maxN)) + (maxn + 0.5 * (n - maxn));
        int tgap = (maxGAP + 0.5 * (GAP - maxGAP)) + 
                   (maxgap + 0.5 * (gap - maxgap));

        int qual = 0;
        int max  = 0;  
        max = MAX(tA,MAX(tC,MAX(tG, MAX(tT, MAX(tN, tgap)))));
        
        if(max == 0){
            uint maxc;
            maxc = MAX(nA,MAX(nC,MAX(nG, MAX(nT, MAX(nN, ngap)))));
            qual = 0;
            if(maxc == nA){
                consensus[i].base = 'A';
            }else if(maxc == nC){
                consensus[i].base = 'C';
            }else if(maxc == nG){
                consensus[i].base = 'G';
            }else if(maxc == nT){
                consensus[i].base = 'T';
            }else if(maxc == nN){
                consensus[i].base = 'N';
            }else if(maxc == ngap){
                consensus[i].base = '-';
            }else{
                fatalf("error in the infer_consenus routine when q = 0");
            }
        }else if(max == tA){   
            consensus[i].base = 'A';
            qual = tA - tC - tG - tT - tN - tgap;
        }else if(max == tC){
            consensus[i].base = 'C';
            qual = tC - tA - tG - tT - tN - tgap;
        }else if(max == tG){
            consensus[i].base = 'G';
            qual = tG - tA - tC - tT - tN - tgap;
        }else if(max == tT){
            consensus[i].base = 'T';
            qual = tT - tA - tC - tG - tN - tgap;
        }else if(max == tN){
            consensus[i].base = 'N';
            qual = tN - tA - tC - tG - tT - tgap;
        }else if(max == tgap){
            consensus[i].base = '-';
            qual = tgap - tA - tC - tG - tT - tN;
        }else{
            fatalf("error in the infer_consenus routine");
        }   
        if(qual < 0) qual = 0;

        consensus[i].qual = qual > 254 ? 255 : qual;
    }
}

/* print the details  of this contig. Return the number of non gap bases in this
 * contig  */
static uint print_contig(const contig* const c, 
                         FILE* const fp,
                         const bool includegaps, 
                         int* const numusefulbases,
                         hashtable* const map)
{
    pre(c != NULL);

    char name[1024];
    sprintf(name, "%d", c->c);
    char* refname = must_find_hashtable(map, name, strlen(name)); 

    /* print the fasta sequence for this contig */
    fprintf(fp, ">Contig%d_%s_%d_%d\n", c->index, refname, c->s, c->e);
    
    uint numcolumns = slcount(c->columns);
    element* consensus = ckalloc(numcolumns * sizeof(element));
    infer_consensus(c, &consensus);

    uint i, j;
    column* col;

    for(i = 0, j = 0, col = c->columns; i < numcolumns; i++, col = col->next){
        if(includegaps == FALSE && consensus[i].base == '-') continue;

        j++;
        if(col->numelements >= 4) *numusefulbases = *numusefulbases + 1;

        fprintf(fp, "%c",consensus[i].base == '-' ? '*' : consensus[i].base);
        if(j % 60 == 0){
            printf("\n");
        }
    }

    if(j % 60 != 0) fprintf(fp, "\n");
    ckfree(consensus);
    return j;
}

/* mark the contigs that are "BAD" i.e. the ones I should not bother printing
 * out to the user */
void mark_bad_contigs(assembly* const assmbl, const int minavgcov)
{
    if(0 == minavgcov) return;

    contig* iter;
    for(iter = assmbl->contigs; iter; iter = iter->next){
        int numbases = 0;
        int numpositions = 0;

        column* c;
        for(c = iter->columns; c; c = c->next){
            numpositions++;
            numbases += c->numelements;
        }

        if(((int)(numbases*1.0/numpositions)) < minavgcov){
            iter->badcontig = TRUE;
        }
    }
}



/* print the details of this assembly */
void print_assembly(const assembly* const a, 
                    FILE* const fp, 
                    const bool includegaps,
                    hashtable* const map)
{    
    pre(a != NULL);

    contig* iter;
    uint numbasescontig UNUSED; 
    
    int numusefulbases = 0;

    uint i;
    for(i = 0, iter = a->contigs; iter; i++, iter = iter->next){
        if(iter->badcontig == TRUE) continue;
        numbasescontig = print_contig(iter,fp,includegaps,&numusefulbases,map);
    }
}

/* print the assembly in the SAM format */
void print_assembly_sam(const assembly* const assmbl,
                        FILE* const samfile,
                        hashtable* const map)
{
    pre(assmbl != NULL);
    pre(samfile != NULL);

    FILE* fp = samfile;

    contig* c;      // a contig in the assembly
    seqread* r;     // a read in a contig
    element* e;     // an element in a read
    
    int rlen;       // length of a read
    char* seq =NULL;// the sequence of a read
    char* qual=NULL;// the quality values in a read

    uint16_t flag;  // for the bitwise flag for each read
    uint pos,start; // 1-based leftmost mapping of the first matching base
    int mq;         // mapping quality, right now we just assign it 255, which
                    // implies that it is not available

    char* rnext = "*";
    int   pnext = 0;
    int   tlen  = 0;

    element* consensus = NULL; // consensus sequence for a contig

    /* print the alignment section now */
    for(c = assmbl->contigs; c; c = c->next){
        if(c->badcontig == TRUE) continue;

        uint numcolumns = slcount(c->columns);

        char name[1024];
        sprintf(name, "%d", c->c);
        char* refname = must_find_hashtable(map, name, strlen(name)); 

        /* the consensus sequence for this contig */
        consensus = ckrealloc(consensus, numcolumns * sizeof(element));
        infer_consensus(c, &consensus);

        uint i = 0, j = 0, k = 0;

        for(i = 0; i < c->numreads; i++){
            r = c->reads[i];
            
            rlen = slcount(r->lunused)+slcount(r->clear)+slcount(r->runused);
            seq  = ckrealloc( seq, (rlen + 1)*sizeof(char));
            qual = ckrealloc(qual, (rlen + 1)*sizeof(char));
            seq[rlen]  = 0;
            qual[rlen] = 0;

            j = 0;
            for(e = r->lunused; e; e = e->next){
                seq[j]  = toupper(e->base);
                qual[j] = e->qual + 33;
                j++;
            }
            for(e = r->clear; e; e = e->next){
                if(e->base == '-') continue;
                seq[j]  = toupper(e->base);
                qual[j] = e->qual + 33;
                j++;
            }
            for(e = r->runused; e; e = e->next){
                seq[j]  = toupper(e->base);
                qual[j] = e->qual + 33;
                j++;
            }
            seq[j] = 0;
            qual[j] = 0;

            flag = 0x0000;
            if(r->complemented == TRUE){
                flag += 0x0010;
            }

            pos   = sldistance(c->columns, r->start);
            forceassert(pos <= numcolumns);
            start = pos;
            for(j = 0, k = 0; j < pos; j++){
                if(consensus[j].base == '-') k++;
            }
            pos -= k;
            mq   = 255;
            
            /* simple way to calculate the mapping quality. Every mismatch from
             * the consensus sequence is punished. */
            for(j = 0, e = r->clear; e; j++, e = e->next){
                forceassert((start + j - 1) < numcolumns);
                if(toupper(consensus[start + j - 1].base) != toupper(e->base)){
                    mq -= e->qual;
                    if(mq < 0) mq = 0;
                }
            }

            fprintf(fp, "%s\t%d\tContig%d_%s_%d_%d\t%d\t%d\t",
                        r->name, flag, c->index, refname, c->s, c->e, pos, mq);

            if(r->lunused != NULL) fprintf(fp, "%dS", slcount(r->lunused));
                
            e = r->clear;
            j = sldistance(c->columns, r->start) - 1;
            while(e){
                forceassert(j < numcolumns);
                /* a match or a mismatch */
                if(e->base != '-' && consensus[j].base != '-'){
                    k = 0;
                    while(e && e->base != '-' && consensus[j].base != '-'){
                        k++;
                        j++;
                        e = e->next;
                    }     
                    fprintf(fp, "%dM", k);
                }else if(e->base == '-' && consensus[j].base != '-'){
                    k = 0;
                    while(e && e->base == '-' && consensus[j].base != '-'){
                        k++;
                        j++;
                        e = e->next;
                    }     
                    fprintf(fp, "%dD", k);
                }else if(e->base != '-' && consensus[j].base == '-'){
                    k = 0;
                    while(e && e->base != '-' && consensus[j].base == '-'){
                        k++;
                        j++;
                        e = e->next;
                    }     
                    fprintf(fp, "%dI", k);
                }else if(e->base == '-' && consensus[j].base == '-'){
                    k = 0;
                    while(e && e->base == '-' && consensus[j].base == '-'){
                        k++;
                        j++;
                        e = e->next;
                    }     
                    fprintf(fp, "%dP", k);
                }
            }            

            if(r->runused != NULL) fprintf(fp, "%dS", slcount(r->runused));

            fprintf(fp, "\t%s\t%d\t%d\t%s\t%s\n", 
                        rnext, pnext, tlen, seq, qual);
        }
    }

    ckfree(consensus);
    ckfree(seq);
    ckfree(qual);
}

/* print the assembly in the ACE format*/
void print_assembly_ace(const assembly* const assmbl, 
                        FILE* const acefile)
{
    pre(assmbl != NULL);
    pre(acefile != NULL);

    FILE* fp = acefile;

    contig* c;               /* details of the contigs */
    int numbases;            /* total number of bases in a contig (incl. *)*/
    seqread* r;              /* a single read */
    element* e;              /* a single element (base) of a read */
    element* consensus = NULL;

    /* the number of elements in a read */
    uint lunused, clear, runused;
       
    int i, j, start;

    /* the first line has the AS segment */ 
    uint total_reads = 0;
    for(c = assmbl->contigs; c; c = c->next){
        total_reads += c->numreads;
    }

    for(c = assmbl->contigs; c; c = c->next){
        if(c->badcontig == TRUE) continue;

        numbases = slcount(c->columns);

        consensus = ckrealloc(consensus, numbases * sizeof(element));
        infer_consensus(c, &consensus);        

        /* the contig line followed by the consensus sequence */
        fprintf(fp, "CO %d %d %d 1 U\n", c->index, numbases, c->numreads);
        for(i = 0; i < numbases; i++){
            fprintf(fp, "%c", consensus[i].base=='-' ? '*' : consensus[i].base);
            if((i + 1) % 60 == 0){
                fprintf(fp, "\n");
            }
        }
        fprintf(fp, "\n\n");

        /* the quality values of the consensus */
        fprintf(fp, "BQ\n");
        for(i = 0; i < numbases; i++){
            if(consensus[i].base != '-'){
                fprintf(fp, "%d ", consensus[i].qual);
            }
            if((i + 1) % 60 == 0){
                fprintf(fp, "\n");
            }
        }
        fprintf(fp, "\n\n");

        /* the AF section, which tells us whether a read is complemented and the
         * offset of the read in the contig */
        for(i = 0; i < (int)c->numreads; i++){
            r = c->reads[i];
            start = sldistance(c->columns, r->start) - slcount(r->lunused);

            fprintf(fp, "AF %s %c %d\n", 
                    r->name, r->complemented ? 'C' : 'U', start);
        }
        fprintf(fp, "\n");

        for(i = 0; i < (int)c->numreads; i++){
            r = c->reads[i];
            
            /* how many bases do I have in total for this read */
            lunused = slcount(r->lunused);
            clear   = slcount(r->clear);
            runused = slcount(r->runused);     
 
            start = lunused + clear + runused;
            fprintf(fp, "RD %s %d 0 0\n", r->name, start);
            
            j = 0;
            for(e = r->lunused; e; e = e->next){
                fprintf(fp, "%c", tolower(e->base));
                if((j + 1) % 60 == 0){
                    fprintf(fp, "\n");
                }
                j++;
            }
            for(e = r->clear; e; e = e->next){
                fprintf(fp, "%c", e->base == '-' ? '*' : toupper(e->base));
                if((j + 1) % 60 == 0){
                    fprintf(fp, "\n");
                }
                j++;
            }
            for(e = r->runused; e; e = e->next){
                fprintf(fp, "%c", tolower(e->base));
                if((j + 1) % 60 == 0){
                    fprintf(fp, "\n");
                }
                j++;
            }
            fprintf(fp, "\n");
            fprintf(fp, "\n");

            fprintf(fp, "QA 1 %d %d %d\n", 
            lunused + clear + runused, lunused + 1 , lunused + clear);
            fprintf(fp, "\n");
        }

        
    }          

    if(consensus != NULL) ckfree(consensus);
}

/* free the resources from this assembly */
void free_assembly(assembly** pa)
{
    assembly* a = *pa;

    contig* iter = a->contigs;
    for(; iter; iter = iter->next){
        column* c = iter->columns;
        for(; c; c = c->next){
            ckfree(c->elements);
        }
        slfreelist(&iter->columns);

        uint i = 0;
        for(; i < iter->numreads; i++){
            seqread* r = iter->reads[i];
            ckfree(r->name);
            slfreelist(&r->lunused);
            slfreelist(&r->clear);
            slfreelist(&r->runused);
            ckfree(r);
        }

        ckfree(iter->reads);    
    }
    slfreelist(&a->contigs);
        
    ckfree(a);
    *pa = NULL;
}

/* calculate the alignment score for this contig. Iterate through the columns
 * and for each column find the consensus. Add 1 for every base in the column
 * that does not agree with the consensus */
static int calculate_alignment_score(const contig* const c)
{  
    pre(c != NULL);

    uint numcolumns = slcount(c->columns);

    element* consensus = ckalloc(numcolumns * sizeof(element));
    infer_consensus(c, &consensus);

    column* col;    
    uint i, j;
    int score = 0;
    for(j = 0, col = c->columns; col; j++, col = col->next){
        for(i = 0; i < col->numelements; i++){
            if(toupper(col->elements[i]->base) != consensus[j].base) score++;
        }
    }

    ckfree(consensus);

    return score;
}

/* remove this read  from the contig. */
static void remove_read(contig* const c, seqread* const r)
{
    pre(c != NULL); 
    pre(r != NULL);
    pre(r->contig == c);

    uint i;
    column* prev;
    column* curr;
    element* p;
    element* e;

    /* remove the read from the contig */
    for(curr = r->start, e = r->clear; 
        curr && curr != r->end->next; 
        curr = curr->next, e = e->next){
        for(i = 0; i < curr->numelements; i++){
            if(curr->elements[i] == e){ 
                curr->elements[i] = NULL;
                break;
            }
        }
    }

    for(p = NULL, e = r->clear; e; p = e, e = e->next){
        if(e && e->base == '-'){
            if(p != NULL){
                p->next = e->next;
                if(e->next != NULL) e->next->prev = p;
            }else{
                r->clear = e->next;
                e->prev = NULL;
            }
            ckfree(e);
            e = p;
        }
    }
    p->next = NULL;

    /* does this change the contig? If there is a column which has all '-', then
     * that should be removed*/
    prev = r->start->prev;
    curr = r->start;
    bool flag;                       // should we remove this column?
    while(curr != NULL && curr->next != NULL){
        flag = TRUE;

        for(i = 0; i < curr->numelements; i++){
            if(curr->elements[i] != NULL && curr->elements[i]->base != '-'){
                flag = FALSE;
                break;
            }
        }  

        if(flag == TRUE){
            /* remove this column */
            if(curr == r->start){
                r->start = prev == NULL ? c->columns : prev;
            }
            if(curr == r->end){
                r->end = curr->next == NULL ? prev : curr->next;
            }

            for(i = 0; i < curr->numelements; i++){
                if(curr->elements[i] != NULL){
                    forceassert(curr->elements[i]->base == '-');
                    curr->elements[i]->prev->next = curr->elements[i]->next;
                    if(curr->elements[i]->next != NULL){
                        curr->elements[i]->next->prev = curr->elements[i]->prev;
                    }
                    ckfree(curr->elements[i]);
                }
            }

            if(prev != NULL) {
                prev->next = curr->next;
                curr->next->prev = prev;
            }else{
                /* this is the first column in the contig */
                forceassert(r->start == c->columns);
                c->columns = curr->next;
                r->start   = curr->next;
                curr->next->prev = NULL; 
            }

            ckfree(curr->elements);
            ckfree(curr);
            
            curr = prev;
        }
        if(curr == r->end) break;

        prev = curr;
        curr = curr == NULL ? r->start : curr->next;
    }
}

static float scoredissimilarity(const column* const el1,const char el2)
{
    pre(el1 != NULL);

    float score = 0;

    uint A = 0, C = 0, G = 0, T = 0, N = 0, GAP = 0;

    uint i = 0;
    for(; i < el1->numelements; i++){
        if(el1->elements[i] != NULL){
            switch(el1->elements[i]->base){
                case 'A': A += 1; break;
                case 'a': A += 1; break;
                case 'C': C += 1; break;
                case 'c': C += 1; break;
                case 'G': G += 1; break;
                case 'g': G += 1; break;
                case 'T': T += 1; break;
                case 't': T += 1; break;
                case 'N': N += 1; break;
                case 'n': N += 1; break;
                case '-': GAP += 1; break;
                default : fatalf("unknown symbol: %c", el1->elements[i]);
            }   
        }   
    }

    char consensus = ' ';
    uint count = MAX(A,MAX(C,MAX(G,MAX(T,MAX(N,GAP)))));
    if(count == A){
        consensus = 'A';
        score += ((C + G + T + N + GAP) * WEIGHTSCORE2 / i);
    }else if(count == C){
        consensus = 'C';
        score += ((A + G + T + N + GAP) * WEIGHTSCORE2 / i);
    }else if(count == G){
        consensus = 'G';
        score += ((A + C + T + N + GAP) * WEIGHTSCORE2 / i);
    }else if(count == T){
        consensus = 'T';
        score += ((A + C + G + N + GAP) * WEIGHTSCORE2 / i);
    }else if(count == N){
        consensus = 'N';
        score += ((A + C + G + T + GAP) * WEIGHTSCORE2 / i);
    }else if(count == GAP){
        consensus = '-';
        score += ((A + C + G + T + N) * WEIGHTSCORE2 / i);
    }else{
        fatal("unknown consensus");
    }
    forceassert(consensus != ' ');

    score += toupper(el2) == consensus ? 0 : WEIGHTSCORE1;

    return score;
}

/* realign this read to the contig */
static void realign_read(contig* const c, seqread* const r)
{
    pre(c != NULL);
    pre(r != NULL);
    forceassert(r->contig == c);

//    fprintf(stdout, "Before aligning %s\n", r->name);
//    print_contig_columns(c);    

    // counters in loops
    uint i, j, k;

    // the read to be realigned; save things for easier access
    uint l1 = slcount(r->clear) + 1; 
    element** Elements = ckallocz(l1 * sizeof(element*));
    element* e;
    for(e = r->clear, i = 1; e; e = e->next, i++){
        Elements[i] = e;
    }
    forceassert(i == l1);

    // the columns to be aligned
    column* rstart = r->start;
    forceassert(rstart != NULL);
    column* rend = r->end;
    forceassert(rend != NULL);

    uint l2 = sldistance(rstart, rend) + 1;
    column** Columns   = ckallocz(l2 * sizeof(column*));
    column* col;
    for(col = rstart, i=1; col && col != rend->next; col = col->next,i++){
        Columns[i] = col;    
    }
    forceassert(i == l2);

    V = ckrealloc(V, (l1 * l2) * sizeof(float));
    for(i = 0; i < l1; i++) V[i] = i;

    I = ckrealloc(I, (l1 * l2) * sizeof(int));

    F = ckrealloc(F, l1 * sizeof(float));
    for(i = 0; i < l1; i++) F[i] = i;

    float G, E;
 
    for(i = 1; i < l2; i ++){
        V[i*l1] = i;

        for(j = 1; j < l1; j++){
            G = scoredissimilarity(Columns[i], Elements[j]->base);
            G += V[(i - 1) * l1 + j - 1];

            E = 1;
            E += V[i * l1 + j - 1];

            F[j] = scoredissimilarity(Columns[i], '-');
            F[j] += V[(i - 1) * l1 + j];

            if(G <= MIN(E,F[j])){
                V[i*l1 + j] = G;
                I[i*l1 + j] = 0;
            }else if(E < F[j]){
                V[i*l1 + j] = E;
                I[i*l1 + j] = 1;
            }else{
                V[i*l1 + j] = F[j];
                I[i*l1 + j] = 2;
            }
            
        }
    }

//    for(i = 0; i < l2; i ++){
//        for(j = 0; j < l1; j++){
//            fprintf(stderr, "%2.2f\t", V[i*l1+j]);
//        }
//        fprintf(stderr, "\n");
//    }
// 
//    for(i = 0; i < l2; i ++){
//        for(j = 0; j < l1; j++){
//            fprintf(stderr, "%d\t", I[i*l1+j]);
//        }
//        fprintf(stderr, "\n");
//    }

    float score UNUSED = V[l1 * l2 - 1];
    i = l2 - 1;
    j = l1 - 1;

    int dir  = 0;
    int flag = -1;

    column* el1;    
    element* el2;
    while((dir =  I[(i * l1) + j]) == 2){
        // remove the spots from this column
        k = 0;
        el1 = Columns[i];
        while(k < el1->numelements){
            if(el1->elements[k] == NULL){
                memmove(el1->elements + k,
                        el1->elements + k + 1,
                       (el1->numelements - k - 1) * sizeof(element*));
                el1->numelements--;
                el1->elements = ckrealloc(el1->elements,
                                          el1->numelements * sizeof(element*));
                k--;
            }
            k++;
        }

        I[(i * l1) + j] = flag;
        i--;
        flag = dir;
        score = V[(i * l1) + j];
    }

    r->end = Columns[i];

    int index;
    forceassert(j == l1 - 1);
    while(j >= 1 && i >= 1){
        dir = I[(i * l1) + j];
        I[(i * l1) + j] = flag;

        if(0 == dir){
            /* the element j-1 aligns with the column i-1 */
            el1 = Columns[i];

            /* do we have an empty spot? */
            index = -1;
            for(k = 0; k < el1->numelements; k++){
                if(el1->elements[k] == NULL){
                    index = k;
                    break;
                }
            }

            if(index == -1){
                el1->numelements++;
                el1->elements = ckrealloc(el1->elements, 
                                        el1->numelements * sizeof(element*));
                el1->elements[el1->numelements - 1] = Elements[j];
            }else{
                el1->elements[index] = Elements[j];
            }

            i--; j--;
        }else if(1 == dir){
            column* prev = Columns[i];

            /* how many of the reads in the previous column do not end there */
            int l = 0;
            for(k = 0; k < prev->numelements; k++){
                if(prev->elements[k] != NULL && 
                   prev->elements[k]->next != NULL) l++;
            }

            /* a new column for the contig */
            column* col = ckalloc(sizeof(column));
            col->numelements = l + 1;
            col->elements = ckalloc(col->numelements * sizeof(element*));
            for(k = 0, l = 0; k < prev->numelements; k++){
                if(prev->elements[k] == NULL || 
                   prev->elements[k]->next == NULL) continue;
                element* tmp = prev->elements[k]->next;

                el2 = ckalloc(sizeof(element));
                el2->base = '-';
                el2->qual = MIN(prev->elements[k]->qual,tmp->qual);

                prev->elements[k]->next = el2;
                el2->prev = prev->elements[k];

                el2->next = tmp;
                tmp->prev = el2;
                
                col->elements[l] = el2;
                l++;
            }
            
            el2 = Elements[j];
            col->elements[l] = el2;

            /* add this column in its location */
            column* tmp = prev->next; 

            if(tmp == NULL || prev == r->end) r->end = col;
            prev->next = col;
            col->prev  = prev;
            col->next  = tmp;
            if(tmp != NULL) tmp->prev = col;

            j--;
        }else if(2 == dir){
            /* a gap in the read at column i-1 */
            el1 = Columns[i];

            index = -1;
            for(k = 0; k < el1->numelements; k++){
                if(el1->elements[k] == NULL){
                    index = k;
                    break;
                }
            }

            if(index == -1){
                el1->numelements++;
                el1->elements = ckrealloc(el1->elements, 
                                          el1->numelements * sizeof(element*));
                index = el1->numelements - 1;
            }

            element* tmp = Elements[j]->next;

            /* the new gap element */
            el2 = ckalloc(sizeof(element));
            el2->base = '-';
            if(tmp == NULL){
                el2->qual =  Elements[j]->qual;
            }else{
                el2->qual = MIN(tmp->qual, Elements[j]->qual);
            }

            /* insert it in place (in the read) */
            Elements[j]->next = el2;
            el2->prev = Elements[j];
            
            el2->next = tmp;
            if(tmp != NULL) tmp->prev = el2;

            /* add it to the column */
            el1->elements[index] = el2;

            i--;
        }else{
            fatalf("incorrect operation in the dp matrix");
        }

        flag  = dir;
        score = V[(i * l1) + j];
    }

    r->start = Columns[i + 1];

    if(i <= 1 && j >= 1){
        // part of the read did not align. Add the non-aligned part to the left
        // over read
        e = sllast(r->lunused);
    
        if(e != NULL){
            e->next = r->clear;
        }else{
            r->lunused = r->clear;
        }
         
        Elements[j+1]->prev->next = NULL;
        Elements[j+1]->prev = NULL;
        r->clear = Elements[j+1];

    }else if(j <= 1 && i >= 1){
        // part of the columns did not align. Just remove the spots in the
        // columns and then change the start for the read.
        r->start = Columns[i + 1];
        for(; i >= 1; i--){
            // remove the spots from this column
            k = 0;
            el1 = Columns[i];
            while(k < el1->numelements){
                if(el1->elements[k] == NULL){
                    memmove(el1->elements + k,
                            el1->elements + k + 1,
                           (el1->numelements - k - 1) * sizeof(element*));
                    el1->numelements--;
                    el1->elements = ckrealloc(el1->elements,
                                    el1->numelements * sizeof(element*));
                    k--;
                }
                k++;
            }
        }
    }

    forceassert(r->lunused == NULL || r->lunused->prev == NULL);
    forceassert(r->clear->prev == NULL);
    forceassert(r->runused == NULL || r->runused->prev == NULL);
    forceassert(r->start != NULL);
    forceassert(r->end != NULL);
    assert(slcount(r->clear) == sldistance(r->start, r->end));

    ckfree(Elements);
    ckfree(Columns);

//    fprintf(stdout, "After aligning %s\n", r->name);
//    print_contig_columns(c);
//    fprintf(stdout, ">>>>>>>>>>>>>>>>> EOR <<<<<<<<<<<<<<<<<<<<<\n");

}

/* realign this contig. */
static void realign_contig(contig* const c)
{
    pre(c != NULL);

    uint i;
    seqread* r;
    column* col;
    bool flag;

    /* no need to realign anything if there is only one read in the contig */
    if(c->numreads == 1) return;

    for(i = 0; i < c->numreads; i++){
        r = c->reads[i];  

        flag = TRUE;
        for(col = r->start; col && col != r->end->next; col = col->next){
            if(col->numelements == 1){
                flag = FALSE;
                break;
            }
        }
        
        if(flag == FALSE) continue;

        /* remove this read from the contig and modify it (remove any columns
         * where we have all bases as gaps)*/ 
        remove_read(c, r);

        /* realign the read to the contig */
        realign_read(c, r);    
    }
    
    uint j, k, numnulls;
    for(col = c->columns; col; col = col->next){
        forceassert(col->numelements >= 1);
        for(i = 0, numnulls = 0; i < col->numelements; i++){
            if(col->elements[i] == NULL){
                numnulls++;
            }
        }

        if(numnulls > 0){
            for(j = 0, k = 0; j < col->numelements; j++){
                if(col->elements[j] != NULL){
                    col->elements[k] = col->elements[j];
                    k++;
                }
            }
            col->numelements = k;
            col->elements = ckrealloc(col->elements,
                                      col->numelements * sizeof(element*));
        }
    }
}

/* realign the assembly */
void realign_assembly(assembly* const assmbl, const bool realignonce)
{
    pre(assmbl != NULL);

    contig* c;    // a single contig in the assembly
    int score;    // score of the contig 
    int newscore; // the score after realignment

    fprintf(stderr, "--------------------------------------------\n");
    for(c = assmbl->contigs; c; c = c->next){
        fprintf(stderr, "Contig%d\n", c->index);

        /* realign this contig until we do not see any improvements in the
         * scores */
        newscore = calculate_alignment_score(c);

        do {
            score = newscore;
            realign_contig(c);
            
            /* what is the new alignment score for this contig*/
            newscore = calculate_alignment_score(c);
            fprintf(stderr, "Old score: %d, New Score: %d\n", score, newscore);
            
            /* do we want to realign just once? */
            if(realignonce == TRUE) break;
        }while(newscore < score);
        fprintf(stderr, "--------------------------------------------\n");
    }    
}

/* free all the static variables in this file. Primarily I dont want to see any
 * more "still reachable" records in valgrind */
void free_alignment_resources()
{   
    ckfree(V);
    ckfree(I);
    ckfree(F);
    ckfree(Columns);
    ckfree(Elements);
}

