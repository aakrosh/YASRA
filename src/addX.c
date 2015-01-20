// addX -- put three Xs at each end of the contigs;
#include "util.h"
#include "seq.h"

int main(int argc, char **argv) {
	SEQ *sf;
	int i;
	
	for (i = 1; i < argc; ++i) {
		sf = seq_open(argv[i]);
		while (seq_read(sf))
			printf("%s\nXXX%sXXX\n", SEQ_HEAD(sf), SEQ_CHARS(sf));
		(void)seq_close(sf);
	}
	return EXIT_SUCCESS;
}
