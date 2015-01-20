/* This module takes in a hits fasta file, which has the alignment information
 * in the header e.g.:

chrM 672 772 HWI-ST407_110127_0082_A80L25ABXX:5:23:2425:180195#0/2 - 0 100 AAAAAAAAAAAAAAAAAAAAAAAAAA  ##################!!!!$$$$$$$$$$^^^^^^^
chrM 13553 13631 HWI-ST407_110127_0082_A80L25ABXX:5:46:7598:187779#0/2 + 0 78 CCCCCCCCCCCCCCCCCCCCC  #################AAAAAAAARRRRRRRHHHHHHH

  The columns are 
    reference name 
    start1
    end1
    read name
    read strand
    start2
    end2
    complete read sequence
    complete read qualities

  YASRA expects just one reference name in the hits file and complains if that
  is not true.

   The output is the consensus assembled sequence. Other optional outputs 
   include ACE format, SAM format and so on
*/

#include "utilities.h"
#include "contig.h"
#include "graph.h"

/* do we want to print the debug information */
int debug_flag = 0;

#if DEBUG
static void check_assembly(const assembly* const assmbl)
{
    pre(assmbl != NULL);

    contig* con;
    column* col;
    for(con = assmbl->contigs; con; con = con->next){
        for(col = con->columns; col; col = col->next){
            assert(col->next == NULL || col->next->prev == col);
        }
    }
}
#endif

static void process_graph(graph* const g, 
                          const bool dorealign,
                          const bool dorealignonce,
                          FILE* const acefile,
                          FILE* const samfile,
                          hashtable* const map,
                          int* const numcontigs,
                          int* const numbases,
                          const int minavgcov)
{
    pre(g != NULL);

    /* we want the reads ordered in the same order as the reference */
    slreverse(&g->node_list);

    /* now we have a graph with the reads as nodes and an edge between two nodes
     * if the reads align. Now lets find a path in the graph which would find
     * the maximum spanning tree. Sort the edges based on their alignment to 
     * the reference genome. The edges with higher weights come first in the 
     * order.*/
    int edgeindex, numedges;
    edge** path = find_maximum_spanning_tree(g, &numedges);
    assert(allnodes_visited(g) == TRUE);
    if(1 == debug_flag){
        printf("Spanning Tree:\n");
        for(edgeindex = 0; edgeindex < numedges; edgeindex++){
            printf("%s -- %s\n", path[edgeindex]->n1->name, 
                                 path[edgeindex]->n2->name);
        }
    }

    /* reads that do not align to anything can be removed right here */
    node* iter = g->node_list;
    for(; iter; iter = iter->next){
        if(0 == iter->num_edges){
            seqread* r = iter->val;
            ckfree(r->name);
            slfreelist(&r->lunused);
            slfreelist(&r->clear);
            slfreelist(&r->runused);
            ckfree(r);
        }
    } 

    /* this is the new assembly and a base with the symbol of a GAP */
    assembly* assmbl = new_assembly();

    /* align the reads in the given order of these edges where the most similar
     * nodes are aligned first */   
    seqread* r1;
    seqread* r2;
    int weight;
    column* start = NULL;

    for(edgeindex = 0; edgeindex < numedges; edgeindex++){  
        /* a point to remember would be that since we use Prims algorithm, one
         * of the nodes from this edge is being added the first time. This can
         * be used to simplify the alignment code. */
        r1 = path[edgeindex]->n1->val;
        r2 = path[edgeindex]->n2->val;
        weight = path[edgeindex]->weight;
        forceassert(r1->c == r2->c);
        start = NULL;

        if(r2->contig == NULL && r1->contig != NULL){
            /* estimate the rough starting of the alignment before you call the
             * routine to actually align it */
            if(r1->s > r2->s){
                uint i;
                for(i = 0; i < r1->contig->numreads; i++){
                    if(r1->contig->reads[i]->s <= r2->s && 
                       r1->contig->reads[i]->e >= r2->s){
                        start = r1->contig->reads[i]->start;
                        break;
                    }
                }
            }else{
                start = r1->start;
            }

            if(start == NULL){
                fprintf(stderr, 
                "Warning: Could result in multiple contigs for the same region\n");
            }else{
                align_nodes(assmbl, r1, start, r2, weight);
            }
        }else if(r1->contig == NULL && r2->contig != NULL){
            if(r2->s > r1->s){
                uint i;
                for(i = 0; i < r2->contig->numreads; i++){
                    if(r2->contig->reads[i]->s <= r1->s &&
                       r2->contig->reads[i]->e >= r1->s){
                        start = r2->contig->reads[i]->start;
                        break;
                    }
                }
            }else{
                start = r2->start;
            }

            if(start == NULL){
                // this can happen if r1 had a hit on the reference where r1->s
                // < r2->s. In that case you would expect r1 to have seeded the
                // contig, rather than r2. 
                fprintf(stderr, 
                "Warning: Could result in multiple contigs for the same region\n");
            }else{
                align_nodes(assmbl, r2, start, r1, weight);
            }
        }else if(r1->contig == NULL && r2->contig == NULL){
            if(r1->s <= r2->s){
                align_nodes(assmbl, r1, NULL, r2, weight);
            }else{
                align_nodes(assmbl, r2, NULL, r1, weight);
            }
        }else{
            fatal("Unhandled condition: both reads have assigned contigs");
        }
#if DEBUG
        check_assembly(assmbl);
#endif
    }

    /* get the contigs in the order of the reference */
    slreverse(&assmbl->contigs);

    /* realign around indels to get better alignments */
    if(dorealign == TRUE){
        realign_assembly(assmbl, dorealignonce);
    }
 
    /* do we print all the contigs we generated in this component? */
    mark_bad_contigs(assmbl, minavgcov);

    /* print the assembled contigs */
    print_assembly(assmbl, stdout, FALSE, map);

    /* print the assembly in the ACE format */
    if(acefile != NULL){
        print_assembly_ace(assmbl, acefile);
    }
    
    /* print the assembly in SAM format */
    if(samfile != NULL){
        print_assembly_sam(assmbl, samfile, map);
    }

    contig* c;
    for(c = assmbl->contigs; c; c = c->next){
        if(c->badcontig == TRUE) continue;
        *numcontigs += 1;
        *numbases   += slcount(c->columns);
    }

    free_assembly(&assmbl);
    ckfree(path);
}

static void process_hits_file(const char* const hitsname,
                              const bool dorealign,
                              const bool dorealignonce,
                              const bool doorient,
                              const bool declone,
                              const bool doconvert,
                              FILE* const acefile,
                              FILE* const samfile,
                              int*  const numcontigs,
                              int*  const numbases,
                              const int minavgcov) 
{
    /* to read the short name from the header */
    char name[1024];
    uint c1 = 0, s1, e1, oldc1 = 0, olds1 = 0, olde1 = 0;
    uint s2, e2;
    char* ptr;

    char* ptr1;  // pointer to the sequence of the read
    char* ptr2;
    char* ptr3;  // pointer to the quality of the read

    /* does the read map reverse complemented */
    bool rc = FALSE;
    char strand;

    /* use these as markers to separate non-contained and contained reads */
    uint end = 0;

    /* variable to read the hitsfile */
    size_t fsize = 1;
    char* fptr = ckalloc(fsize + 1);
    FILE* fp = ckopen(hitsname, "r");
    
    int count = 0;  // counter for the reads

    /* the graph structure that we intend to utilize */
    graph* g = new_graph();
    node*  n = NULL;

    /* a map to put the index and the name of the reference sequence */
    hashtable* map = new_hashtable(10);
    char indexname[1024];

    /* print the header in the samfile */
    if(samfile != NULL) fprintf(samfile, "@HD\tVN:1.3\n");

    char* chromname = NULL; // the reference sequence
    char refname[1024];     // the name of the reference sequence for one read

    while(getline(&fptr, &fsize, fp) != -1){
        if(*fptr == '#') continue;

        if(sscanf(fptr, "%s\t%d\t%d\t%[^'\t']\t%c\t%d\t%d\t%*s\t%*s\n", 
                        refname, &s1, &e1, name, &strand, &s2, &e2) == 7){

            if(chromname == NULL){
                chromname = copy_string(refname);

                sprintf(indexname, "%d", c1);
                forceassert(
                lookup_hashtable(map, indexname, strlen(indexname)) == NULL);
                add_hashtable(map, 
                              indexname, strlen(indexname), 
                              copy_string(chromname));
            }
            
            if(strncmp(chromname, refname, strlen(chromname)) != 0){
                ckfree(chromname);
                chromname = copy_string(refname);
                c1 += 1;

                sprintf(indexname, "%d", c1);
                forceassert(
                lookup_hashtable(map, indexname, strlen(indexname)) == NULL);
                add_hashtable(map, 
                              indexname, strlen(indexname), 
                              copy_string(chromname));
            }                  

            // set the strand information
            rc  = FALSE;
            if(strand == '-') rc = TRUE;
        }

        // ignore this read if this is a clone of the last read
        if(declone == TRUE && 
           c1 == oldc1 && s1 == olds1 && e1 == olde1) continue;

        /* get the sequence and quality information (if that is available) */
        int slen = 0;
        int i = 0;
        ptr2  = NULL;
        ptr   = fptr;

        while(*ptr != '\n'){
            if(*ptr == '\t'){
                i++;
                while(*ptr == '\t' && *ptr != '\n')ptr++;
            }
            if(i == 7) break;
            ptr += 1;
        }
        ptr1 = ptr;

        while(*ptr != '\n'){
            if(*ptr == '\t'){
                ptr2 = ptr;
                i++;
                while(*ptr == '\t' && *ptr != '\n')ptr++;
            }
            if(i == 8) break;
            ptr += 1;
        }
        ptr3 = ptr;
        
        if(*ptr3 == '\n'){
            // no quality values this time
            ptr2 = ptr3;
            ptr3 = NULL;
        }

        // the length of the sequence 
        forceassert(ptr1 != NULL);
        forceassert(ptr2 != NULL);
        slen = ptr2 - ptr1;
        forceassert(slen > 0);

        n = new_node(name,
                     new_read(ptr1,
                              ptr3,
                              slen,
                              name,
                              rc,
                              c1,
                              s1, e1,
                              s2, e2,
                              doorient,
                              doconvert));

        
        /* is this new node a part of the current contig. If we know this will
         * be the beginning of a new contig, then we should probably stop here
         * and process the current graph. */
        if((c1 != oldc1) || (end != 0 && s1 > end)){
            process_graph(g, dorealign, dorealignonce, acefile, samfile, map, numcontigs, numbases, minavgcov);
            free_graph(&g);

            /* prepare for the next set of reads */
            g = new_graph();
            if(c1 != oldc1){
                olds1 = 0;
                olde1 = 0;
                end = 0;
            }
        }

        if(g->node_list != NULL) forceassert(c1 == oldc1);
        end = e1 > end ? e1 : end;

        /* add the node to the current graph */
        add_node(g, n);    

        if(++count % 1000 == 0){
            fprintf(stderr, "Added %d reads to the graph\n", count);
        }        
        
        /* lets remember the locations where this read aligns to the template */
        oldc1 = c1;
        olds1 = s1;
        olde1 = e1;
    }

    if(g->node_hash->elcount > 0){
        process_graph(g, dorealign, dorealignonce, acefile, samfile, map, numcontigs, numbases, minavgcov);
    }

    free_hashtable_completely(&map);

    free_graph(&g);
    ckfree(fptr);
    if(chromname != NULL) ckfree(chromname);
    fclose(fp);
}

static void assembler(const char* const hitsname,
                      const char* const acename,
                      const char* const samname,
                      const bool dorealign,
                      const bool dorealignonce,
                      const bool doorient,
                      const bool declone,
                      const bool doconvert,
                      const int minavgcov)
{
    /* lets print the input details to the assembler */
    fprintf(stderr, "Hits file: %s\n", hitsname);
    if(acename != NULL)    fprintf(stderr, "Output ACE file: %s\n", acename);
    if(samname != NULL)    fprintf(stderr, "Output SAM file: %s\n", samname);
    if(dorealign == TRUE)  fprintf(stderr, "Will realign the reads\n");     
    fprintf(stderr, "--------------------------------------------\n");

    FILE* acefile = NULL;
    if(acename != NULL){
        acefile = ckopen(acename, "w");
        fprintf(acefile, "AS ");
        int i = 0;
        for(; i < 100; i++) fprintf(acefile, " ");
        fprintf(acefile, "\n");
    }

    FILE* samfile = NULL;
    if(samname != NULL){
        samfile = ckopen(samname, "w");
    }

    /* Iterate through the hits file. Once you have found the reads that
     * seem like they would be part of the same contig, create a graph for those
     * reads, find the maximum spanning tree. Align the reads in order as per
     * the maximum spanning tree and then write the details in the output files
     * . Then repeat the process for the next contig */
    int numcontigs = 0;
    int numbases   = 0;

    process_hits_file(hitsname, 
                      dorealign, 
                      dorealignonce,
                      doorient, 
                      declone, 
                      doconvert, 
                      acefile, 
                      samfile,
                      &numcontigs,
                      &numbases,
                      minavgcov);

    if(acename != NULL){
        /* prepend the AS segment to the file. Some viewers like Tablet do not *
         * work in the absence of this segment */
        if(fseek(acefile, 3, SEEK_SET) != 0){
            fatal("error in prepending header to the ACE file");
        }

        fprintf(acefile, "%d %d", numcontigs, numbases);            

        fclose(acefile);
    }


    if(samname != NULL) fclose(samfile);

    timestamp("Finished assembly. ");
}
 
static void usage()
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, 
"\t-d,--debug: if set, print debug information (slow)\n");
    fprintf(stderr, 
"\t-r,--realign: if set, realign the reads (slow)\n");
    fprintf(stderr, 
"\t-i,--illumina: if set, quality values are in Illumina format [sanger]\n"); 
    fprintf(stderr, 
"\t-o,--orient: if set, nucs2 in the LASTZ output needs to be oriented\n");
    fprintf(stderr,
"\t-c,--declone: if set, remove clones before assembling\n");
    fprintf(stderr,
"\t-h,--hits: the hits file with the LASTZ output\n");
    fprintf(stderr,
"\t-a,--ace: the name of the output ACE file\n");
    fprintf(stderr,
"\t-s,--sam: the name of the output SAM file\n");
    return;
}

int main(int argc, char** argv)
{
    argv0 = "assembler";
    int c;

    /* the input files?*/
    char* hitsname     = NULL;
    char* acename      = NULL;
    char* samname      = NULL;

    /* options */
    bool dorealign     = FALSE;
    bool doorient      = FALSE; 
    bool declone       = FALSE;
    bool doconvert     = FALSE;
    bool dorealignonce = FALSE;

    /* ignore contigs with average coverage lower than this */
    int minavgcov = 0;

    while (1){
        static struct option long_options[] = {
            {"debug"   , no_argument      , 0, 'd'},
            {"hits"    , required_argument, 0, 'h'},
            {"ace"     , required_argument, 0, 'a'},
            {"sam"     , required_argument, 0, 's'},
            {"cov"     , required_argument, 0, 'm'},
            {"realign" , no_argument      , 0, 'r'},
            {"realign1", no_argument      , 0, '1'},
            {"orient"  , no_argument      , 0, 'o'}, 
            {"declone" , no_argument      , 0, 'c'},
            {"illumina", no_argument      , 0, 'i'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc,argv,"dr1iocm:h:a:s:",long_options,&option_index);

        if (c == -1) break;

        switch (c){
            case 0:
                break;
            case 'd':
                debug_flag = 1;
                break;
            case 'h':
                hitsname  = optarg;
                break;
            case 'a':
                acename = optarg;
                break;
            case 's':
                samname = optarg;
                break;
            case 'r':
                dorealign = TRUE;
                break;
            case '1':
                dorealign     = TRUE;
                dorealignonce = TRUE;
                break;
            case 'o':
                doorient = TRUE;
                break;
            case 'c':
                declone = TRUE;
                break;
            case 'i':
                doconvert = TRUE;
                break;
            case 'm':
                minavgcov = atoi(optarg);
            case '?':
                break;
            default:
                abort();
        }
    }

    forceassert(optind == argc);

    /* we need the hits file */
    if(hitsname == NULL){
        fatal("error : please provide a hits file as input");
        usage();
    }

    /* allocate resources for heavy recursion in this module */
    allocate_resources();

    /* start clock book-keeping */
    t0 = time(0);

    /*make sure the files are there */
    struct stat buf;
    if(stat(hitsname, &buf) != 0){
        fatalf("error in finding the hits file:%s (%s)",hitsname,strerror(errno));
        usage();
    }    
    
    assembler(hitsname, 
              acename, 
              samname,
              dorealign, 
              dorealignonce,
              doorient,
              declone,
              doconvert,
              minavgcov);

    /* this is primarily for Valgrind. Relinquish all the static variables in
     * contig.c. */
    free_alignment_resources();

    /* print the relevant stats used by the program */
    print_usage();

    return EXIT_SUCCESS;
}
