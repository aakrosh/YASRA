/* trim_assembly -- remove places where there are <= min hits or at least half
*   as many rejects as hits.
*
*   Syntax:
	trim_assembly assembly.fa hits.fa reject.fa [min=3]
* where assembly.fa is created by make_assembly (to glue contigs together).
*/

#include "util.h"
#include "seq.h"

#define MAX_READ_LEN 10000
#define MIN_CONTIG_LEN 30

int *hit, *reject;
char buf[MAX_READ_LEN];

void mark(char *filename, int *x) {
	FILE *fp = ckopen(filename, "r");
	int i, b, e;

	while (fgets(buf, MAX_READ_LEN, fp))
		if (buf[0] == '>') {
			if (sscanf(buf, "%*s %d %d", &b, &e) != 2)
				fatalf("endpoints: %s", buf);
			for (i = b; i < e; ++i)
				++x[i];
		}
	fclose(fp);
}

int main(int argc, char **argv) {
	SEQ *sf;
	uchar *s;
	int i, b, e, N, c;

	argv0 = "trim_assembly";
	int min_cov = 2;

	while(argc > 4){
		argc--;
		if(strncmp(argv[argc],"min=", 4)){
			min_cov = atoi(argv[argc]+4);	
		}
	}

	if (argc != 4)
		fatal("args: assembly.fa hits.fa reject.fa");
	sf = seq_get(argv[1]);
	s = SEQ_CHARS(sf);
	N = SEQ_LEN(sf);
	hit = ckalloc(N*sizeof(int));
	reject = ckalloc(N*sizeof(int));
	for (i = 0; i < N; ++i)
		hit[i] = reject[i] = 0;
	mark(argv[2], hit);
	mark(argv[3], reject);
	for (i = 0; i < N; ++i)
		if (hit[i] <= min_cov || 2*reject[i] >= hit[i])
			s[i] = 'N';
	for (b = c = e = 0; b < N; b = e) {
		// move b past a run of "N"
		while (b < N && s[b] == 'N')
			++b;
		if (b == N)
			break;
		if (b > e)
			fprintf(stderr, "discarded: %d-%d\n", e, b);
		// move e to the end of this contig
		for (e = b; e < N && s[e] != 'N'; ++e)
			;
		if (e-b > MIN_CONTIG_LEN) {
			fprintf(stderr, "Contig%d %d %d\n", c, b, e);
			printf(">Contig%d\n", c++);
			for (i = b; i < e; ++i)
				putchar(s[i]);
			putchar('\n');
		} else
			fprintf(stderr,
			  "discarded short contig %d-%d\n", b, e);
	}
	return 0;
}
