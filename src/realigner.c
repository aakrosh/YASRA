/* This file contains Eric's test program.
   This version uses a linked list to implement the algorithm exactly,
   culls out columns of -'s after each iteration,
   uses a combined scoring scheme,
   and uses a constant band size around the old path for the pairwise 
     alignment step. And it's hopefully 'optimized to the max.'
   Added support for the 'N' character.
 
   1.01 - Made program work for any size input
   1.02 - Simplified column structure by attaching header for frags
   1.03 - Took out printing of input
   1.04 - Added consensus
   1.05 - bug fix, increased size of Decode[]
   1.06 - allows input of all characters, unusual chars treated as n's
	  strings are also now printed in the same case they are inputed
   1.07 - bug fix, order of adjusting i & j in insert switched
   1.08 - took out using variable to declare array size
   1.09 - Added constant CONF to control captitalizing of consensus
   1.10 - took out last array declared with a variable size, added '.'s to
	  alignment printout
   2.00 - Allow for multiple contigs
   2.01 - Fix on capital for consensus with 'n's, dots in output for single
	  frag
   2.02 - Allow for comments and calculatate % mismatch
   3.00 - Changed i/o to vertical form
   3.01 - added MIX define type and bug fix to i/o
   3.02 - accept multiple blank lines between contigs & fixed possible
	  infinite loop problem.
   3.03 - fixed error when reAligned frag ends with an insert
   3.04 - Changed to take lowest scoring path that ends as close as possible
	  to the original path.
   3.05 - Made the bandsize an optional input
   3.06 - Aakrosh: added LINE macro and corresponding changes. Changed colInf 
      and colDepth to be int. Changed '-' to '_'.
   3.07 - Aakrosh: made changes to accomodate ACE files.
*/


#define MIX .5
#define LINE 1024000

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "util.h"

/* run `cat MAlign | wc -L` to fix this */
int	MaxDepth = 1024000;
/* run `cat hits.fa | grep '>' | wc -l` to fix this*/
int	MaxFragNum = 512000;
/* run `cat hits.fa | grep '>' | wc -L` to fix this */
int	MaxFragSize = 5120;

int	NumIts = 0;
int	NumCons = 0;

int	DoneFlag = 0;

int	BandSize = 8;

int	Rows;

FILE* ace = NULL;
FILE* tmpFile = NULL;
char* coInfo = NULL;

typedef struct fragEl {
    char		el;
    struct fragEl	*prev;
    struct fragEl *next;
    struct fragEl *up;
    struct fragEl *down;
}
fragEl;

typedef struct colStruc {
    int			colInf[6];
    int			colDepth;
    double		preCalc[6];
    struct colStruc	*prev;
    struct colStruc	*next;
    fragEl		frags;
}
colStruc;

typedef struct {
    fragEl	*ofrag;
    fragEl	*lastEl;
    char		*efrag;
    char		*afrag;
    int		maxLen;
    int		len;
    colStruc	*scol;
    colStruc	*ecol;
    int		row;
	/*these members are for the ace files*/
	char* name;				/*name of the read*/
	char complement;		/*'C' if it is complemented, 'U' otherwise*/
	char* seq;				/*the whole sequence as is, doesnt change*/
	char* qa;				/*the QA line,doesnt change*/
}
frag;

frag		*Frags;
colStruc	*FirstCol;
colStruc	*LastCol;

fragEl		*FreeElList=NULL;

char *morespace(char *arr, int *top, int factor) {
    char *t;
    int   prev;

    t = arr;
    prev = *top;
    *top = prev + (prev >> 1) + 8;
    arr = malloc(factor * (*top));
    if (prev != 0) {
        arr = memcpy(arr, t, prev * factor);
        free(t);
    }
    return(arr);
}

fragEl *getEl() {
    fragEl	*ptr;
    register int	i;

    if (FreeElList == NULL) {
        FreeElList = ptr = (fragEl *) malloc(sizeof(fragEl)*4096);
        for (i = 0; i < 4095; ++i) {
            ptr->next = ptr + 1;
            ++ptr;
        }
        ptr->next = NULL;
    }

    ptr = FreeElList;
    FreeElList = ptr->next;
    return ptr;
}

void freeEl(ptr) fragEl *ptr;
{
    ptr->next = FreeElList;
    FreeElList = ptr;
}

colStruc *FreeColList=NULL;

colStruc *getCol() {
    register int	i;
    colStruc	*ptr;

    if (FreeColList == NULL) {
        FreeColList = ptr = (colStruc *) malloc(sizeof(colStruc)*4096);
        for (i = 0; i < 4095; ++i) {
            ptr->next = ptr+1;
            ++ptr;
        }
        ptr->next = NULL;
    }

    ptr = FreeColList;
    FreeColList = ptr->next;
    for (i=0; i < 6; ++i)
        ptr->colInf[i] = 0;
    ptr->colDepth = 0;
    ptr->preCalc[5] = 0.0;
    ptr->next = ptr->prev = NULL;
    ptr->frags.up = ptr->frags.down = &(ptr->frags);
    ptr->frags.prev = (fragEl *) &ptr;
    return ptr;
}

void freeCol(ptr) colStruc *ptr;
{
    ptr->next = FreeColList;
    FreeColList = ptr;
}


int		encode[128];

int readFrags() {
    register char *s;
    register int	i, r, numFrags;
    char		buffer[LINE];
    int		*row, temp;
    colStruc	*curCol;
    fragEl	*elPtr, *telPtr;

    row = (int *) malloc(MaxDepth*sizeof(int));
    buffer[LINE-1] = '\0';
    numFrags = 0;
    for (i=0; i < MaxDepth; ++i)
        row[i] = -1;
    elPtr = telPtr = getEl();
    elPtr->el = -1;
    for (i=1; i < BandSize; ++i) {
        elPtr->next = getEl();
        elPtr = elPtr->next;
        elPtr->el = -1;
    }

    FirstCol = curCol = getCol();

    Rows = -1;

	unsigned long n = 1;
	char* lineptr = ckalloc(n+1);
	char name[MAX_SEQ_NAME];

	if(ace != NULL &&
	   getline(&lineptr, &n, ace) != -1){
	   	assert(strncmp(lineptr,"CO", 2) == 0);
		coInfo = copy_string(lineptr);
	}

    while ((s=fgets(buffer, LINE-2, stdin)) != NULL) {
        i = strlen(buffer);
        if (buffer[i-1] == '\n')
            buffer[--i] = '\0';
        if (i > LINE-2) {
            fprintf(stderr,"Each input line must not be more than %d chars\n", LINE);
            exit(1);
        }
        if (buffer[0] == '%') {
            printf("%s\n",buffer);
            continue;
        }
        r = 0;
        curCol->next = getCol();
        curCol->next->prev = curCol;
        curCol = curCol->next;
        for ( ; *s != '\0'; ++s) {
            if (*s != ' ') {
                if (r > Rows)
                    Rows = r;
                if (r >= MaxDepth) {
                    row = (int *) morespace((char *) row, &MaxDepth, sizeof(int));
                    for (i=Rows; i < MaxDepth; ++i)
                        row[i] = -1;
                }
                if (row[r] == -1) {
					/*we have a new read*/
                    row[r] = numFrags;
					
					/*lets read the info from the ace file*/
					if(ace != NULL && getline(&lineptr, &n, ace) != -1){
						assert(strncmp(lineptr, "AF", 2) == 0);

						if(sscanf(lineptr, "AF %s %c 1", name, &Frags[numFrags].complement) == 2){
							Frags[numFrags].name = copy_string(name);
							if(getline(&lineptr, &n, ace) != -1){
								assert(strncmp(lineptr,"RD", 2) == 0);
							}
							if(getline(&lineptr, &n, ace) != -1){
								Frags[numFrags].seq = copy_string(lineptr);
							}
							if(getline(&lineptr, &n, ace) != -1){
								assert(strncmp(lineptr,"QA", 2) == 0);
								Frags[numFrags].qa = copy_string(lineptr);
							}
						}
					}

                    Frags[numFrags].scol = curCol;
                    Frags[numFrags].row = r;
                    Frags[numFrags].maxLen = MaxFragSize;
                    Frags[numFrags].efrag =
                        (char *) malloc(MaxFragSize*sizeof(char));
                    Frags[numFrags].afrag =
                        (char *) malloc(MaxFragSize*sizeof(char));
                    Frags[numFrags].len = 0;
                    elPtr = getEl();
                    for (i=0; i < BandSize; ++i) {
                        elPtr->el = -1;
                        elPtr->next = getEl();
                        elPtr->next->prev = elPtr;
                        elPtr = elPtr->next;
                    }
                    Frags[numFrags].ofrag = Frags[numFrags].lastEl = elPtr;
                    ++numFrags;
                    if (numFrags == MaxFragNum)
                        Frags = (frag *) morespace((char *)Frags,&MaxFragNum,
                                                   sizeof(frag));
                } else
                    elPtr = Frags[row[r]].lastEl;
                if ((i = encode[(int)*s]) < 0) {
                    fprintf(stderr,"Illegal char in input line %d\n",r);
                    exit(0);
                }
                ++(curCol->colInf[i]);
                ++(curCol->colDepth);
                elPtr->el = i;
                elPtr->up = &(curCol->frags);
                elPtr->down = curCol->frags.down;
                curCol->frags.down = elPtr;
                elPtr->down->up = elPtr;
                elPtr->next = Frags[row[r]].lastEl = getEl();
                Frags[row[r]].lastEl->prev = elPtr;
                if (i != 0) {
                    Frags[row[r]].afrag[Frags[row[r]].len] = *s;
                    Frags[row[r]].efrag[Frags[row[r]].len++] = i;
                    if (Frags[row[r]].len == Frags[row[r]].maxLen) {
                        temp = Frags[row[r]].maxLen;
                        Frags[row[r]].afrag = morespace(Frags[row[r]].afrag,
                                                        &temp,sizeof(char));
                        Frags[row[r]].efrag = morespace(Frags[row[r]].efrag,
                                                        &Frags[row[r]].maxLen,sizeof(char));
                        if (MaxFragSize < Frags[row[r]].maxLen)
                            MaxFragSize = Frags[row[r]].maxLen;
                    }
                }
            } else
                if (row[r] != -1) {
                    Frags[row[r]].lastEl->el = -1;
                    Frags[row[r]].lastEl->next = telPtr;
                    Frags[row[r]].lastEl = Frags[row[r]].lastEl->prev;
                    Frags[row[r]].ecol = curCol->prev;
                    row[r] = -1;
                }

            ++r;
        }
        while (r <= Rows) {
            if (row[r] != -1) {
                Frags[row[r]].lastEl->el = -1;
                Frags[row[r]].lastEl->next = telPtr;
                Frags[row[r]].lastEl = Frags[row[r]].lastEl->prev;
                Frags[row[r]].ecol = curCol->prev;
                row[r] = -1;
            }
            ++r;
        }
        if (curCol->colDepth == 0)
            break;
    }
    if (s == NULL)
        DoneFlag = 1;
    free(row);
    curCol->next = LastCol = getCol();
    LastCol->prev = curCol;
    ++Rows;
	ckfree(lineptr);
    return numFrags;
}


int		Decode[6];


void printAlign(numFrags) int numFrags;
{

    int		row, i;
    fragEl	**curEl;
    char		**curChar;
    colStruc	*col, *ecol;
    fragEl	*ptr, *kep;
    void		*tmp;

    col = FirstCol;
    while (col != LastCol)
        if (col->colDepth == col->colInf[0]) {
            if (col != FirstCol)
                col->prev->next = col->next;
            else
                FirstCol = col->next;
            col->next->prev = col->prev;
            ptr=col->frags.down;
            for (i=0; i < col->colDepth; ++i) {
                ptr->prev->next = ptr->next;
                ptr->next->prev = ptr->prev;
                kep = ptr;
                ptr = ptr->down;
                freeEl(kep);
            }
            ecol = col;
            col = ecol->next;
            freeCol(ecol);
        } else {
            col->frags.el = 0;
            col = col->next;
        }
    for (i=0; i < numFrags; ++i)
        Frags[i].scol->frags.el = 1;
    curEl = (fragEl **) malloc(sizeof(fragEl *)*Rows);
    curChar = (char **) malloc(sizeof(char *)*Rows);

	frag** pFrags = (frag **) ckallocz(sizeof(frag *)* Rows);

	if(tmpFile){
		fprintf(tmpFile,"%s", coInfo);
	}

    for (i=0; i < Rows; ++i) {
        curEl[i] = NULL;
        curChar[i] = NULL;
    }
    col = FirstCol;
    while (col != LastCol) {
        if (col->frags.el)
            for (i=0; i < numFrags; ++i)
                if (col == Frags[i].scol) {
                    if (curEl[Frags[i].row] == NULL) {
						pFrags[Frags[i].row] = &Frags[i];
                        curEl[Frags[i].row] = Frags[i].ofrag;
                        curChar[Frags[i].row] = Frags[i].afrag;
                    } else {
                        tmp = (void *) curEl;
                        curEl = (fragEl **) malloc(sizeof(fragEl *)*(Rows+1));
                        memcpy(curEl, tmp, sizeof(fragEl *) * Rows);
                        free(tmp);
                        tmp = (void *) curChar;
                        curChar = (char **) malloc(sizeof(char *)*(Rows+1));
                        memcpy(curChar, tmp, sizeof(char *) * Rows);
                        free(tmp);
                        tmp = (void *) pFrags;
						pFrags = (frag **) malloc(sizeof(frag **)*(Rows+1));
                        memcpy(pFrags, tmp, sizeof(char *) * Rows);
						free(tmp);
						curEl[Rows] = Frags[i].ofrag;
                        curChar[Rows] = Frags[i].afrag;
						pFrags[Rows] = &Frags[i];
                        ++Rows;
                    }
                }
		if(tmpFile){		
			for(i =0; i < Rows; i++){
				if(pFrags[i] != NULL){
					/*write the ace file*/
					fprintf(tmpFile, "AF %s %c %d\n", pFrags[i]->name, pFrags[i]->complement, i+1);
					fprintf(tmpFile, "RD %s L 0 0\n", pFrags[i]->name);
					fprintf(tmpFile, "%s", pFrags[i]->seq);
					fprintf(tmpFile, "%s", pFrags[i]->qa);
					pFrags[i] = NULL;
				}
			}
		}


        for (row=0; row < Rows; ++row) {
            if (curEl[row] == NULL)
                printf(" ");
            else if (curEl[row]->el == -1) {
                printf(" ");
                curEl[row] = NULL;
                curChar[row] = NULL;
            } else {
                if (curEl[row]->el == 0)
                    printf("_");
                else
                    printf("%c",*curChar[row]++);
                curEl[row] = curEl[row]->next;
            }
        }
        printf("\n");
        col = col->next;
    }

    printf("\n");
    free(curEl);
    free(curChar);
    return;
}

void printModel() {
    int		i;
    colStruc	*j;

    printf("\n");
    for (j = FirstCol; j != LastCol; j = j->next) {
        for (i = 0; i < 6; ++i)
            printf("%d   ",j->colInf[i]);
        printf("\n");
    }
    for (i = 0; i < 6; ++i)
        printf("%d   ",LastCol->colInf[i]);
    printf("\n");
}

#define DEL 1
#define SUB 0
#define INS 2


double	*mat;
char	*bmat;
int	*bmatPtr;
int	BmatSize;
int	*shift;


/* prcedure reAlign realigns fragment number fnum against the alignment of
   thre rest of the fragments */

void reAlign (fnum) int fnum;
{
    register int	i, j, m, n;
    char		*fel, *cptr;
    double	min, dval, sval, ival;
    double	*lval, *tval;
    frag		*pf;
    int		mlen, max, mark;
    colStruc	*col, *minCol, *tcol, *mstart, *mstop;
    fragEl	*ptr, *tptr;
    fragEl	*fptr;

    pf = &Frags[fnum];
    mark = 0;


    /* Strip fragment from structure */
    col = pf->scol;
    for (ptr=pf->ofrag; ptr->el != -1; ptr=ptr->next) {
        ptr->up->down = ptr->down;
        ptr->down->up = ptr->up;
        --(col->colDepth);
        --(col->colInf[(int)ptr->el]);
        col = col->next;
    }

    mstart = pf->scol;
    mlen = 1+2*BandSize;
    for (n=0; n < BandSize; ++n) {
        if (mstart == FirstCol) {
            FirstCol = mstart->prev = getCol();
            FirstCol->next = mstart;
        }
        mstart = mstart->prev;
    }
    for (j=0, tptr=ptr=pf->ofrag; j < BandSize; ++j) {
        ptr = ptr->next;
        while (ptr->el == 0) {
            ++mlen;
            ptr = ptr->next;
        }
        tptr = tptr->prev;
    }
    for (j=1; j <= mlen; ++j)
        mat[j] = 0.0;

    for (mstop=pf->ecol->next,j=0; mstop!=LastCol && j<BandSize;
            ++j,mstop=mstop->next) ;
    for (col = mstart; col != mstop; col = col->next) {
        m = col->colDepth - col->colInf[5];
        if (m != 0) {
            max = col->colInf[0];
            for (i=1; i < 5; ++i)
                if (max < col->colInf[i])
                    max = col->colInf[i];
            min = m;
            for (i=0; i < 5; ++i) {
                col->preCalc[i] = MIX*(1.0-(double)col->colInf[i]/min);
                if (col->colInf[i] != max)
                    col->preCalc[i] += (1-MIX);
            }
        } else {
            for (i=0; i < 5; ++i)
                col->preCalc[i] = 0.0;
        }
    }

    fel = pf->efrag;
    for (i = 1; i <= pf->len; ++i, ++fel) {
        ptr = ptr->next;
        while (ptr->el == 0) {
            ++mlen;
            mat[mlen] = MaxFragSize;
            ptr = ptr->next;
        }
        mat[mlen+1] = MaxFragSize;
        shift[i] = 1;
        while (tptr->el == 0) {
            --mlen;
            ++(shift[i]);
            tptr = tptr->next;
            mstart = mstart->next;
        }
        tptr = tptr->next;
        col = mstart;
        mstart = mstart->next;
        bmatPtr[i] = mark;
        cptr = &bmat[mark];
        mark += mlen;
        if (mark > BmatSize) {
            bmat = morespace(bmat,&BmatSize,sizeof(double));
            cptr = &bmat[mark-mlen];
        }
        tval = mat;
        lval = &mat[shift[i]];
        for (j=1; j <= mlen && col != LastCol; ++j, col=col->next) {
            dval = (*tval++) + col->preCalc[0];
            sval = (*lval++) + col->preCalc[(int)*fel];
            if (*fel != 5)
                ival = *lval+1.0;
            else
                ival = *lval;
            if (sval <= dval && sval <= ival) {
                mat[j] = sval;
                *cptr = SUB;
            } else if (dval <= ival) {
                mat[j] = dval;
                *cptr = DEL;
            } else {
                mat[j] = ival;
                *cptr = INS;
            }
            ++cptr;
        }
    }

    cptr = &bmat[bmatPtr[pf->len]];
    for (n=1,col=mstart->prev;col!=pf->ecol && n<=mlen; ++n,col=col->next) ;
    min = mat[n];
    minCol = col;
    cptr = &bmat[bmatPtr[pf->len]];
    j = n+1;
    col = minCol->next;
    tcol = minCol->prev;
    for (i=n-1; i > 0; --i) {
        if (j <= mlen && col != LastCol) {
            if (mat[j] < min || (mat[j] == min && cptr[n-1]==DEL)) {
                n = j;
                min = mat[j];
                minCol = col;
            }
            ++j;
            col = col->next;
        }
        if (mat[i] < min || (mat[i] == min && cptr[n-1]==DEL)) {
            n = i;
            min = mat[i];
            minCol = tcol;
        }
        tcol = tcol->prev;
    }

    ptr = pf->lastEl;
    mlen = j = n-1;
    i = pf->len;
    fel = &(pf->efrag[pf->len-1]);
    col = minCol;

    while  (i > 0) {
        if (bmat[bmatPtr[i]+j] == SUB) {
            ptr->el = m = *fel--;
            ++(col->colDepth);
            ++(col->colInf[m]);
            ptr->up = &(col->frags);
            ptr->down = col->frags.down;
            ptr->up->down = ptr;
            ptr->down->up = ptr;
            col = col->prev;
            j = j+shift[i]-1;
            --i;
        } else if (bmat[bmatPtr[i]+j] == DEL) {
            ptr->el = 0;
            ++(col->colDepth);
            ++(col->colInf[0]);
            ptr->up = &(col->frags);
            ptr->down = col->frags.down;
            ptr->up->down = ptr;
            ptr->down->up = ptr;
            col = col->prev;
            --j;
        } else {
            tcol = getCol();
            tcol->prev = col;
            tcol->next = col->next;
            col->next->prev = tcol;
            col->next = tcol;
            ++(tcol->colDepth);
            ptr->el = m = *fel--;
            ++(tcol->colInf[m]);
            ptr->down = ptr->up = &(tcol->frags);
            tcol->frags.down = tcol->frags.up = ptr;
            tcol->frags.prev = (fragEl *) tcol;
            fptr = col->frags.down;
            for (n=0; n < col->colDepth; ++n) {
                if (fptr->next->el != -1) {
                    ++(tcol->colDepth);
                    ++(tcol->colInf[0]);
                    tptr = getEl();
                    tptr->prev = fptr;
                    tptr->next = fptr->next;
                    tptr->next->prev = tptr;
                    tptr->prev->next = tptr;
                    tptr->el = 0;
                    tptr->up = &(tcol->frags);
                    tptr->down = tcol->frags.down;
                    tptr->up->down = tptr;
                    tptr->down->up = tptr;
                }
                fptr = fptr->down;
            }
            j = j+shift[i];
            --i;
        }
        if (ptr == pf->ofrag && i > 0) {
            pf->ofrag = getEl();
            pf->ofrag->prev = ptr->prev;
            ptr->prev->next = pf->ofrag;
            ptr->prev = pf->ofrag;
            pf->ofrag->next = ptr;
        }
        ptr = ptr->prev;
    }
    pf->ofrag = ptr->next;
    while (ptr->el != -1) {
        ptr->el = -1;
        ptr = ptr->prev;
    }
    if (col != NULL)
        pf->scol = col->next;
    else
        pf->scol = FirstCol;
    if (bmat[bmatPtr[pf->len]+mlen] == INS)
        pf->ecol = (colStruc *) pf->lastEl->down->prev;
    else
        pf->ecol = minCol;

    return;
}


void useErr(char *name) {
    fprintf(stderr,"usage: %s [-b#]\n",name);
    exit(1);
}


int main(argc, argv) int argc;
char *argv[];
{
    register int	i;
    int		numFrags, flag;
    int		score, max, oldScore, n;
    colStruc	*col;
    fragEl	*ep, *kep;
	char* acefile = NULL;

	while(argc > 1){
		--argc;
		if(0 == strncmp(argv[argc],"-ace=", 5)){
			acefile = copy_string(argv[argc]+5);
			ace = ckopen(acefile, "r");
			if(!(tmpFile = tmpfile()) ){
				fprintf(stderr, "error in creating a temp stream");
				exit(1);
			}
		}else if(0 == strncmp(argv[argc],"-b", 2)){
			BandSize = atoi(argv[argc]+2);
			if( BandSize < 1){
            	fprintf(stderr,"Illegal band size\n");
               	exit(1);
			}
		}else{
			useErr(argv[0]);
		}
	}

    /* initialize values */
    Decode[0] = '_';
    Decode[1] = 'a';
    Decode[2] = 'c';
    Decode[3] = 'g';
    Decode[4] = 't';
    Decode[5] = 'n';
    for (i = 0; i < 128; i++)
        encode[i] = 5;
    encode['_'] = 0;
    encode['a'] = encode['A'] = 1;
    encode['c'] = encode['C'] = 2;
    encode['g'] = encode['G'] = 3;
    encode['t'] = encode['T'] = 4;
    encode['n'] = encode['N'] = 5;

    while (!DoneFlag) {
        Frags = (frag *) malloc(MaxFragNum*sizeof(frag));
        numFrags = readFrags();
        if (numFrags == 0)
            continue;
        fprintf(stderr, "Read %d reads\n", numFrags);
        fflush(stderr);

        /* initialize matrices needed for calculations */
        BmatSize = MaxFragSize*(4*BandSize+2);
        bmat = (char *) malloc(BmatSize*sizeof(char));
        mat = (double *) malloc(MaxFragSize*sizeof(double));
        mat[0] = MaxFragSize;
        bmatPtr = (int *) malloc(MaxFragSize*sizeof(int));
        shift = (int *) malloc(MaxFragSize*sizeof(int));

        /* tack 2 blank columns on end of alignment to allow movement */
        LastCol->next = getCol();
        LastCol->next->prev = LastCol;
        LastCol = LastCol->next;
        LastCol->next = getCol();
        LastCol->next->prev = LastCol;
        LastCol = LastCol->next;

        score = 0;
        for (col=FirstCol; col != LastCol; col=col->next) {
            if (col->colDepth > 0) {
                max = col->colInf[0];
                for (i = 1; i < 5; ++i)
                    if (col->colInf[i] > max)
                        max = col->colInf[i];
                score += (col->colDepth-max);
            }
        }

        oldScore = score+1;
        flag = 0;
        while (oldScore > score) {
            oldScore = score;
            ++flag;
            for (i=0; i < numFrags; i++) {
                reAlign(i);
                /*
                printAlign(numFrags);
                */
            }

            score = 0;
            n = 0;
            for (col=FirstCol; n < BandSize; ++n,col=col->next) {
                if (col->colDepth > 0) {
                    max = col->colInf[0];
                    for (i = 1; i < 5; ++i)
                        if (col->colInf[i] > max)
                            max = col->colInf[i];
                    score += (col->colDepth-max);
                }
            }
            for ( ; col != LastCol; col=col->next) {
                max = col->colInf[0];
                if (col->colDepth == max)  /* if column of blanks, remove */
                {
                    col->prev->next = col->next;
                    col->next->prev = col->prev;
                    ep=col->frags.down;
                    for (i=0; i < col->colDepth; ++i) {
                        ep->prev->next = ep->next;
                        ep->next->prev = ep->prev;
                        kep = ep;
                        ep = ep->down;
                        freeEl(kep);
                    }
                }
                else {
                    for (i = 1; i < 5; ++i)
                        if (col->colInf[i] > max)
                            max = col->colInf[i];
                    score += (col->colDepth-max);
                }
            }
        }

        free(bmat);
        free(mat);
        free(bmatPtr);
        free(shift);

        printAlign(numFrags);
        fprintf(stderr,"After %d iterations\n",flag);
        NumIts += flag;
        ++NumCons;

        for (i=0; i < numFrags; ++i) {
            Frags[i].lastEl->next = FreeElList;
            FreeElList = Frags[i].ofrag;
            free(Frags[i].afrag);
            free(Frags[i].efrag);
        }
        free(Frags);
        LastCol->next = FreeColList;
        FreeColList = FirstCol;
    }
    fprintf(stderr,"Total %d its in %d contigs for ave of %5.2f its/con\n",
            NumIts,NumCons,(double)NumIts/(double)NumCons);

	if(ace){
		fclose(ace);
		ace = ckopen(acefile, "w");
		unsigned long n = 1;
		char* lineptr = ckalloc(n+1);
		rewind(tmpFile);
		while(getline(&lineptr, &n, tmpFile) != -1){
			fprintf(ace, "%s", lineptr);
		}

		fclose(ace);
	}

    return EXIT_SUCCESS;
}


