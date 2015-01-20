/* template_hits -- process blastz output for a "template" sequence
*   vs. a FastA file of "reads" (short sequences).
*
* Produces a fasta file of the trimmed hits, but skips trimming where
* the template has a gap (run of Ns).
*
* syntax :
	blastz template.fa 454.fa | template_hits template.fa 454.fa
*
* Optional command-line arguments:
*	debug - report certain intermediate values
* 	stats - produce a table of differences between template.fa and 454.fa
*	cov=70 - require that a specified percentage of the read be aligned
*	pct=80 - require a minumum percent identity
*	discard=filename - save aligned positions FastA headers of weakly
*	  aligning reads (i.e., failing the cov= or pct= criteria)
*/

#include "util.h"
#include "seq.h"

//ignore gaps less than MIN_GAP
#define MIN_GAP 30

// positions of gaps in the template
#define MAX_GAP 100000
struct gap {
	int b, e;
} G[MAX_GAP];

// hits extended into gaps, including the discarded hits
#define MAX_HIT 30000000
struct hit {
	char *header;	// FastA header
	uchar *s;	// nucleotides that align or extend into a gap
	int b1, b2, e1, e2, // alignment ends: template (b1,e1) and read (b2,e2)
	  nmatch, npair,    // number of matches, aligned pairs
	  rc,	// is the hit to the reverse complement of the read?
	  low,	// if reject != 0, this holds the offending percentage
	  extend,	// extention? 0 = no, 1 = left, 2 = right
	  reject; // 0 = keep, 1 = low coverage, 2 = low identity
} H[MAX_HIT];

uchar *template, *aread;	// sequences being compared
int rc, len, b1, b2, e1, e2, nhit, ngap, debug, offset, pct, cov, template_len,
   nread, cov_bp, tot_len, more_template, more_read, template_bp, all_len,
   temp_npair, temp_nmatch, many[4][4], temp_many[4][4], stats, code[256],
   temp_more_template, temp_more_read;
char header[500], *template_header, *nuc = "ACGT";

// for sorting hits based on starting location in the template
// (or ending positions for hits starting at the same position)
int compar(const void *a, const void *b) {
	int i = ((struct hit *)a)->b1 - ((struct hit *)b)->b1,
	    j = ((struct hit *)a)->e1 - ((struct hit *)b)->e1;

	return (i != 0 ? i : j);
}

// record information about an aligned pair of nucleotides
void aligned_pair(char x, char y) {

	++temp_npair;
	if (x == y)
		++temp_nmatch;
	if (strchr(nuc, x) && strchr(nuc, y))
		temp_many[code[(unsigned)x]][code[(unsigned)y]]++;
}

// process the gap-free segments of a local alignment
void get_stats() {
	char buf[500];
	int b1, b2, e1, e2, i, j, old_e1 = 0,  old_e2 = 0;

	for (i = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			temp_many[i][j] = 0;
	temp_nmatch = temp_npair = temp_more_template = temp_more_read = 0;
	while (fgets(buf, 500, stdin) && buf[0] != '}') {
		if (sscanf(buf, "  l %d %d %d %d", &b1, &b2, &e1, &e2) !=4)
			fatalf("line: %s", buf);
		if (old_e1 > 0 && b1 > old_e1+1)
			temp_more_template += (b1 - old_e1 - 1);
		if (old_e2 > 0 && b2 > old_e2+1)
			temp_more_read += (b2 - old_e2 - 1);
		old_e1 = e1;
		old_e2 = e2;
		// blastz-output uses base-1 sequence coordinates
		for (i = b1-1, j=b2-1; i < e1; ++i, ++j)
			aligned_pair(toupper(template[i]), toupper(aread[j]));
	}
	H[nhit].npair = temp_npair;
	H[nhit].nmatch = temp_nmatch;
}

void print_stats() {
	int i, j, N, M;

	fprintf(stderr,
	  "\n#template bases label rows, read bases label columns\n");
	fprintf(stderr, "#        A       C       G       T\n");
	for (M = N = i = 0; i < 4; ++i) {
		fprintf(stderr, "#%c", nuc[i]);
		M += many[i][i];
		for (j = 0; j < 4; ++j) {
			fprintf(stderr, " %7d", many[i][j]);
			if (i != j)
				N += many[i][j];
		}
		fputc('\n', stderr);
	}
	fprintf(stderr,
	  "#%d matches, %d mismatches, %d more in template, %d more in read\n",
	  M, N, more_template, more_read);

	fprintf(stderr, "#%d non-N aligned nucleotides at %4.2f%% identity\n",
	  M+N, 100.0*(float)M/(float)(M+N));
	fprintf(stderr, "%d total differences\n",
	  N + more_template + more_read);
}

// scan the template; record start and (non-inclusive) end positions of gaps
// ignore gaps of length < MIN_GAP
void get_gaps() {
	int i, j;

	for (i = 0; i < template_len; i = j) {
		if (ngap >= MAX_GAP)
			fatal("Too many gaps in the template.");
		while (i < template_len && template[i] != 'N')
			++i;
		if (i == template_len)
			break;
		for (j = i+1; template[j] == 'N'; ++j)
			;
		if (j >= i + MIN_GAP) {
			G[ngap].b = i;
			G[ngap].e = j;
			fprintf(stderr, "gap %d: %d-%d\n", ngap++, i, j);
		}
	}
}

// Extend a hit into an adjacent gap, if one exists, then record it.
// Also handle the X's added by the addX program, if any.
// Determine if coverage is inadequate.
void hit() {
	int i, j, g;

	if (nhit >= MAX_HIT)
		fatal("too many hits");
	// convert to base-0 addresses
	--b1;
	--b2;
	if (debug) {
		printf("hit b1 = %d, b2 = %d, e1 = %d, e2 = %d\n  %s\n",
		  b1, b2, e1, e2, header);
		printf("  rc = %d, read = %s\n", rc, aread);
	}

	// find the closest gap to the left of where the read is aligned
	for (g = -1; g < ngap-1 && G[g+1].e <= b1; ++g)
		;

	// if the read's left end extends into a gap ..
	if (g >= 0 && G[g].e >= b1-3) {
		// extend match to the left end of the read, omitting Xs
		b1 += (offset-b2);
		b2 = offset;
		H[nhit].extend = 1;
	}

	++g;	// go to the next gap
	// reality check -- gap must be to the right of the hit
	if (g < ngap && G[g].b < e1)
		fatalf("impossible: gap %d-%d, hit %d-%d, gap %d-%d\n  %s",
		  (g > 0 ? G[g-1].b : 0), (g > 0 ? G[g-1].e : 0),
		  b1, e1, G[g].b, G[g].e, header);
	// if the read's right end extends into a gap ..
	if (g < ngap && G[g].b <= e1+3) {
		// extend match to the right end of the read, omitting Xs
		e1 += (len-offset-e2);
		e2 = len-offset;
		H[nhit].extend = 2;
	}

	if (100*(e2-b2) < cov*(len-2*offset)) { // even when extended into gaps
		j = 100*(e2-b2)/(len-2*offset);
		if (debug)
			printf("low coverage %d: %s\n", j, header);
		H[nhit].reject = 1;
		H[nhit].low = j;
	}

	if (debug)
		printf("  adjusted b1 = %d, e1 = %d, b2 = %d, e2 = %d\n",
		  b1, e1, b2, e2);
	H[nhit].header = copy_string(header);
	H[nhit].rc = rc;
	H[nhit].b1 = b1;
	H[nhit].e1 = e1;
	H[nhit].b2 = b2;
	H[nhit].e2 = e2;

	if (!H[nhit].reject) {
		H[nhit].s = ckalloc((e2-b2+1)*sizeof(uchar));
		for (j = 0, i = b2; i < e2; ++j, ++i)
			H[nhit].s[j] = aread[i];
		H[nhit].s[e2-b2] = '\0';
		if (debug)
			printf("  aligned string: %s\n", H[nhit].s);
		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j)
				many[i][j] += temp_many[i][j];
		cov_bp += (e2-b2);
		tot_len += (len - 2*offset);
		more_template += temp_more_template;
		more_read += temp_more_read;
	}
}

int main(int argc, char **argv) {
	SEQ *tf, // template
	    *rf; //reads
	FILE *discard = NULL;
	char buf[500], *p, *q;
	int i, k, low_cov, low_pct, nsave;

	argv0 = "template_hits";
	while (argc > 3) {
		if (same_string(p=argv[argc-1], "debug"))
			debug = 1;
		else if (same_string(p, "stats"))
			stats = 1;
		else if (strncmp(p, "cov=", 4) == 0)
			cov = atoi(p+4);
		else if (strncmp(p, "pct=", 4) == 0)
			pct = atoi(p+4);
		else if (strncmp(p, "discard=", 8) == 0)
			discard = ckopen(p+8, "w");
		else
			fatalf("bad arg: %s", p);
		--argc;
	}
	if (argc != 3)
	   fatal("args: template 454.fa [pct=??] [cov=??] [discard=??] [debug] [stats]");

	// get template
	tf = seq_get(argv[1]);
	template = SEQ_CHARS(tf);
	template_header = SEQ_HEAD(tf);
	template_len = SEQ_LEN(tf);
	// template_bp is the number of non-N, non-X characters in the template
	for (template_bp = i = 0; i < template_len; ++i)
		if (!strchr("XN", template[i]))
			++template_bp;
	get_gaps();

	// look at the first read
	rf = seq_get(argv[2]);
	aread = SEQ_CHARS(rf);
	// a shameless hack to handle possible use of addX,
	// avoids having offset as a commond-line option
	for (offset = 0; offset < SEQ_LEN(rf) && aread[offset] == 'X';
	  ++offset)
		;
	fprintf(stderr, "offset %d, cov %d, pct %d, debug %d, stats %d\n",
	  offset, cov, pct, debug, stats);
	nread = 1;
	all_len = SEQ_LEN(rf);

	// translate nucleotide into row/column index in table of aligned pairs
	code['A'] = 0; code['C'] = 1; code['G'] = 2; code['T'] = 3;

	// read the blastz output file from stdin
	if (fgets(buf, 500, stdin) == NULL || !same_string(buf, "#:lav\n"))
		fatalf("%s is not a blastz output file", argv[2]);
	while (fgets(buf, 500, stdin)) {
		if (same_string(buf, "s {\n")) {
			// record the read's length and aligned orientation
			if (fgets(buf, 500, stdin) == NULL ||
			    fgets(buf, 500, stdin) == NULL ||
			    sscanf(buf, "%*s %*s %d %d", &len, &rc) != 2)
				fatal("cannot find revcomp");
		} else if (same_string(buf, "h {\n")) {
			// get the read's fasta header line
			if (fgets(buf, 500, stdin) == NULL ||
			    fgets(buf, 500, stdin) == NULL ||
			    (p = strchr(buf, '"')) == NULL ||
			    (q = strchr(p+1, '"')) == NULL)
				fatal("cannot parse input");
			*q = '\0';
			strcpy(header, p+1);
			// find the read in the file of all reads
			while (!strstr(header, SEQ_HEAD(rf))) {
				if (!seq_read(rf))
					fatalf("bad EOF in %s", argv[2]);
				++nread;
				all_len += SEQ_LEN(rf);
			}
			if (SEQ_LEN(rf) != len)
				fatal("lengths do not agree");
			if (rc)
				rf = seq_revcomp_inplace(rf);
			aread = SEQ_CHARS(rf);
		} else if (same_string(buf, "a {\n")) {
			// find end-points of the alignments
			if (fgets(buf, 500, stdin) == NULL ||
			    fgets(buf, 500, stdin) == NULL ||
			    sscanf(buf, "  b %d %d", &b1, &b2) != 2 ||
			    fgets(buf, 500, stdin) == NULL ||
			    sscanf(buf, "  e %d %d", &e1, &e2) != 2)
				fatal("could not parse alignment info");
			get_stats();  // extract info from gap-free segments
			if (100*temp_nmatch < pct*temp_npair) {
				k = 100*temp_nmatch/temp_npair;
				if (debug)
					printf("low identity %d: %s\n",
					  k, header);
				H[nhit].reject = 2;
				H[nhit].low = k;
			}
			hit();
			++nhit;
		}
	}

	if (debug)
		printf("nhit = %d\n---------------------------\n", nhit);
	if (stats)
		print_stats();
	// sort the hits by start position in the template
	qsort((void *)H, nhit, sizeof(struct hit), compar);
	// print the sorted hits
	for (nsave = low_cov = low_pct = i = 0; i < nhit; ++i) {
		if (H[i].reject == 0) {
			if (H[i].extend == 1)
				printf(">Left");
			else if (H[i].extend == 2)
				printf(">Right");
			else 
				printf(">Hit");
			printf("%d %d %d %s, pos %d-%d\n%s\n", i,
			   H[i].b1, H[i].e1, H[i].header,
			   H[i].b2-offset, H[i].e2-offset, H[i].s);
			++nsave;
		} else if (discard != NULL) {
			if (H[i].reject == 1)
				fprintf(discard, ">Low_cov%d", low_cov++);
			else {
				if (H[i].reject != 2)
					fatal("impossible H[i].reject");
				fprintf(discard, ">Low_pct%d", low_pct++);
			}
			fprintf(discard, ":%d %d %d %s\n",
			  H[i].low, H[i].b1, H[i].e1, H[i].header);
		}
	}

	while (seq_read(rf)) {
		++nread;
		all_len += SEQ_LEN(rf);
	}

	fprintf(stderr, "#%d of %d (%4.2f%%) reads saved\n",
	  nsave, nread, 100.0*(float)nsave/(float)nread);
	if (!stats) {
		int j, M, N;
		for (M = N = i = 0; i < 4; ++i) {
			M += many[i][i];
			for (j = 0; j < 4; ++j) {
				if (i != j)
					N += many[i][j];
			}
		}
		fprintf(stderr,
		  "#%d non-N aligned nucleotides at %4.2f%% identity\n",
		  M+N, 100.0*(float)M/(float)(M+N));
	}
	fprintf(stderr, "#average length = %4.1f, average aligned = %4.1f",
	  (float)tot_len/(float)nsave, (float)cov_bp/(float)nsave);
	fprintf(stderr, " (all reads: %4.1f)\n", (float)all_len/(float)nread);
	if (discard != NULL)
		fprintf(stderr, "#%d < %d%% covered, %d < %d%% identical\n",
		   low_cov, cov, low_pct, pct);
	fprintf(stderr,
	 "#%d bp in aligned segments = %3.1fX coverage of the %d-bp template\n",
	  cov_bp, (float)cov_bp/(float)template_bp, template_bp);
	return 0;
}
