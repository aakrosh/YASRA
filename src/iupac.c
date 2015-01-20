#include "iupac.h"

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

extern const char iupac_code[16];
const char iupac_code[] = {'Z','T','G','K','C','Y','S','B','A','W','R','D','M','H','V','N'};

/*return the numeric code for the base*/
int return_code(const char base)
{
	int A = 8;
	int C = 4;
	int G = 2;
	int T = 1;
	int val = -1;

	switch(base){
		case 'A': 
			val = A; break;
		case 'C': 
			val = C; break;
		case 'G': 
			val = G; break;
		case 'T': 
			val = T; break;
		case 'M': 
			val = A|C; break;
		case 'R': 
			val = A|G; break;
		case 'W': 
			val = A|T; break;
		case 'S': 
			val = C|G; break;
		case 'Y': 
			val = C|T; break;
		case 'K': 
			val = G|T; break;
		case 'V': 
			val = A|C|G; break;
		case 'H': 
			val = A|C|T; break;
		case 'D': 
			val = A|G|T; break;
		case 'B': 
			val = C|G|T; break;
		case 'N': 
			val = A|C|G|T; break;
		default :
			fatalf("unknown base:%c",base);
	}
	return val;
}

/*return the IUPAC base for the following set*/
char get_alleles(const char base, const bool A, const bool C, const bool G, const bool T)
{
	if(base == 'N' ||
	   base == '_' ||
	   base == '*'){
	   return base;
	}

	int nA = A ? 1 : 0;	
	int nC = C ? 1 : 0;	
	int nG = G ? 1 : 0;	
	int nT = T ? 1 : 0;	

	nA = nA << 3;
	nC = nC << 2;
	nG = nG << 1;

	assert((nA|nC|nG|nT) != 0);
	assert((nA|nC|nG|nT) < 16);
	return iupac_code[nA|nC|nG|nT];
}
