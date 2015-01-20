/*This module has some routines for reading completely formatted ACE file as
 * per the schema described on the url: 
 * http://bozeman.mbt.washington.edu/consed/consed.html#documentation
 *
 * This however is not for reading partially written ACE files by realigner 
 * or the assembler module. Modules that need to read/write partial ACE files
 * should do so outside this module.
 */

#ifndef ACE_H
#define ACE_H

#include "util.h"
#include "slist.h"
#include "hash.h"
#include "contig.h"
#include "iupac.h"

/*This is what a gap looks like in the ACE file*/
#define GAP '*'

/*call the consensus based on the bases in the column*/
char call_consensus(const COLUMN* const column, const bool alleles);

/*read the next contig and fill in the requisite datastructures*/
CONTIG* read_contig(FILE* af, char** plineptr, size_t* pn, char** pconsensus, char** pname, struct hash* const hash);
#endif

