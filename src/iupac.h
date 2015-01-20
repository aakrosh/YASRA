#ifndef IUPAC_H
#define IUPAC_H

#include "util.h"

/*return the IUPAC nucleotide symbol 
 * A -> A
 * C -> C
 * G -> G
 * T -> T
 * A || C -> M
 * A || G -> R
 * A || T -> W
 * C || G -> S
 * C || T -> Y
 * G || T -> K
 * A || C || G -> V
 * A || C || T -> H
 * A || G || T -> D
 * C || G || T -> B
 * A || C || G || T -> N
 */

/*return the numeric code for the base*/
int return_code(const char base);

/*return the IUPAC base for the following set*/
char get_alleles(const char base, const bool A, const bool C, const bool G, const bool T);
#endif
