// find substitution spectrum outside of specified hypervariable region
#include "util.h"
#include "seq.h"

#define MAX_NAME_LEN 1000

int AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT;

int count_same(uchar *s, int i, char c) {
	int n, j;

	for (n = 0, j = i; j > 0 && toupper(s[j]) == c; --j, ++n)
		;
	if (n == 0 || s[i] == c)
		while (toupper(s[++i]) == c)
			++n;
	return n;
}
	

int main(int argc, char **argv) {
	argv0 = "substitutions";
	SEQ *sf, *mf;
	uchar *s = NULL, *m, a, b, x, y;
	
	/*variables used by aakrosh to add name of the contig to the report*/
	char name[MAX_NAME_LEN];		/*header of the contig*/
	char* nm = NULL;				/*pointer to the name*/
	char* ter = NULL;				/*terminating character*/
	/*variables end*/
	
	FILE *fp;
	int old_e1 = 0, old_e2 = 0;
	char old_seq_name[500], buf[500], c;
	int b1, e1, b2, e2, i, j, k, n1, n2, rc, N, p, q, match,
	    mismatch, more_argv3, more_argv2, B, E, gaps = 0, substitutions = 0,
	    offset = 0;

	B = 0;
	E = -1;
	while (argc > 1) {
	  	if (same_string(argv[argc-1], "substitutions"))
			substitutions = 1;
		else if (same_string(argv[argc-1], "gaps"))
			gaps = 1;
		else if (strncmp(argv[argc-1], "offset=", 7) == 0)
			offset = atoi(argv[argc-1]+7);
		else
			break;
		--argc;
	}
	if (argc == 6) {
		B = atoi(argv[4]);
		E = atoi(argv[5]);
		argc = 4;
		fprintf(stderr, "%d-%d is hypervariable\n", B, E);
	}
	if (argc != 4)
		fatal("args: blastz.out template reads");
	old_seq_name[0] = '\0';
	mf = seq_get(argv[2]);
	m = SEQ_CHARS(mf) - 1;
	sf = seq_open(argv[3]);
	fp = ckopen(argv[1], "r");
	match = mismatch = more_argv3 = more_argv2 = 0;
	while (fgets(buf, 500, fp)) {
		if (same_string(buf, "s {\n")) {
			if(fgets(buf, 500, fp) == NULL) fatalf("error in reading");
			if(fgets(buf, 500, fp) == NULL) fatalf("error in reading");
			if (sscanf(buf, "%*s %*s %*s %d", &rc) != 1)
				fatalf("rc: %s", buf);
			old_e1 = old_e2 = 0;
		} else if (same_string(buf, "h {\n")) {
			if (!seq_read(sf))
				fatalf("bad EOF on %s", argv[3]);
			if(fgets(buf, 500, fp) == NULL) fatalf("error in reading");
			if(fgets(buf, 500, fp) == NULL) fatalf("error in reading");
			for (k = 0; k < 4 &&!strstr(buf, SEQ_HEAD(sf)); ++k)
				if (!seq_read(sf))
					fatalf("unexpected EOF on %s", argv[3]);
			if (k == 4) {
				fprintf(stderr, "'%s'\n", buf);
				fprintf(stderr, "'%s'\n", SEQ_HEAD(sf));
				fatalf("cannot find %s", SEQ_HEAD(sf));
			}
			if (rc == 1)
				sf = seq_revcomp_inplace(sf);
			s = SEQ_CHARS(sf) - 1;
			if(strlen(SEQ_HEAD(sf)) > MAX_NAME_LEN){
				fatal("Increase the width of the line");
			}
			strcpy(name, SEQ_HEAD(sf));
			nm = name+1;
			if((ter = strchr(nm,' '))!= NULL){
				*ter = '\0';
			}	
		} else if (same_string(buf, "a {\n"))
			old_e1 = old_e2 = 0;
		else if (sscanf(buf, "  l %d %d %d %d",
			    &b1, &b2, &e1, &e2) == 4) {
			if (old_e2 > 0 && (b1 < B || e1 > E)) {
				p = b1 - old_e1 - 1;
				q = b2 - old_e2 - 1;
				if (p > 0 && q > 0)
			  fatalf("b1 = %d, old_e1 = %d, b2 = %d, old_e2 = %d\n",
			    b1, old_e1, b2, old_e2);
				if (p <= 0 && q <= 0)
					fatal("neither p or q is positive");
				if (p > 0) {
					if (gaps) {
					    c = toupper(m[old_e1+1]);
					    n1 = count_same(m, old_e1+1, c);
					    n2 = count_same(s, old_e2, c);
					    if (n1 - n2 == p) {
						printf("%s:%d %d %c, ",
						  argv[2], old_e1, n1, c);
						printf("%s.%s:%d %d %c\n",
						  argv[3], nm, old_e2-offset,
						  n2, c);
					    } else for (k = 1; k <= p; ++k)
						printf("%s:%d extra %c\n",
						  argv[2], old_e1+k, 
						  toupper(m[old_e1+k]));
					}
					more_argv2 += p;
				} else {
					if (gaps)
					    for (k = 1; k <= q; ++k)
						printf("%s.%s:%d extra %c\n",
						  argv[3], nm, old_e2+k-offset,
						  toupper(s[old_e2+k]));
					more_argv3 += q;
				}
			}
			old_e1 = e1;
			old_e2 = e2;
			for (i = b1, j = b2; i <= e1; ++i, ++j) {
				if (i >= B && i <= E)
					continue;
				a = m[i];
				x = toupper(a);
				b = s[j];
				y = toupper(b);
				if (substitutions && x != y &&
					x != 'N' && y != 'N')
					printf("%s:%d %c -> %c %s.%s:%d\n",
					  argv[2], i, a, b, argv[3], nm, j-offset);
				if (x == 'A' && y == 'A')
					++AA;
				else if (x == 'A' && y == 'C')
					++AC;
				else if (x == 'A' && y == 'G')
					++AG;
				else if (x == 'A' && y == 'T')
					++AT;
				else if (x == 'C' && y == 'A')
					++CA;
				else if (x == 'C' && y == 'C')
					++CC;
				else if (x == 'C' && y == 'G')
					++CG;
				else if (x == 'C' && y == 'T')
					++CT;
				else if (x == 'G' && y == 'A')
					++GA;
				else if (x == 'G' && y == 'C')
					++GC;
				else if (x == 'G' && y == 'G')
					++GG;
				else if (x == 'G' && y == 'T')
					++GT;
				else if (x == 'T' && y == 'A')
					++TA;
				else if (x == 'T' && y == 'C')
					++TC;
				else if (x == 'T' && y == 'G')
					++TG;
				else if (x == 'T' && y == 'T')
					++TT;
			}
		}
	}

	printf("%s bases label rows, %s bases label columns\n",
	  argv[2], argv[3]);
	printf("       A      C      G      T\n");
	printf("A %6d %6d %6d %6d\n", AA, AC, AG, AT);
	printf("C %6d %6d %6d %6d\n", CA, CC, CG, CT);
	printf("G %6d %6d %6d %6d\n", GA, GC, GG, GT);
	printf("T %6d %6d %6d %6d\n", TA, TC, TG, TT);
	printf("%d matches, %d mismatches, %d more in %s, %d more in %s\n",
	  AA + CC + GG + TT,
	  (i = AC + AG + AT + CA + CG + CT + GA + GC + GT + TA + TC + TG),
	  more_argv2, argv[2], more_argv3, argv[3]);

	N = AA+AC+AG+AT+CA+CC+CG+CT+GA+GC+GG+GT+TA+TC+TG+TT;
	printf("%d aligned nucleotides at %4.2f%% identity\n",
	  N, 100.0*(float)(AA + CC + GG + TT)/(float)N);
	printf("%d total differences\n", i + more_argv2 + more_argv3);
	return 0;
}
