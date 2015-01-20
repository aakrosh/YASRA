/* make_template -- separate contigs by a run of Ns to make a template
*
* syntax: make_template contigs.fa [noends] [N=??] [info=filename] [min=?]
*
* noends -- do not put a run of Ns at each end
* N=37   -- insert runs of 37 Ns
* info=filename	-- write lines like "Contig3 12345 21543" to the named file
*    each line shows the location of a contig in the template
* min=40 -- discard contigs with less than 40 bp
*/

#include "util.h"
#include "seq.h"

#define MIN_CONTIG_LEN 150

int Ns=1000, cur_pos;

void put_Ns() {
	int i;

	for (i = 0; i < Ns; ++i)
		putchar('N');
	putchar('\n');
	cur_pos += Ns;
}

int main(int argc, char **argv) {
	SEQ *sf;
	FILE *fp = NULL;
	char *p;
	int k, n, ends = 1, min_len = MIN_CONTIG_LEN;

	argv0 = "make_template";
	while (argc > 2) {
		if (same_string(p = argv[argc-1], "noends"))
			ends = 0;
		else if (strncmp(p, "N=", 2) == 0)
			Ns = atoi(p+2);
		else if (strncmp(p, "info=", 5) == 0)
			fp = ckopen(p+5, "w");
		else if (strncmp(p, "min=", 4) == 0)
			min_len = atoi(p+4);
		else
			fatalf("improper command-line argument: %s", p);
		--argc;
	}
	if (argc != 2)
		fatal("arg: contigs.fa [noends] [N=?] [info=filename] [min=?]");
	sf = seq_open(argv[1]);
	printf("> make_template\n");
	for (n = 0; seq_read(sf); n = 1)
		if ((k = SEQ_LEN(sf)) >= min_len) {
			if (ends || n > 0)
				put_Ns();
			if (fp != NULL)
				fprintf(fp, "Contig%d %d %d\n",
				  n, cur_pos, cur_pos + k);
			cur_pos += k;
			printf("%s\n", SEQ_CHARS(sf));
		}
	if (ends)
		put_Ns();
	if (fp != NULL)
		fclose(fp);
	return 0;
}
