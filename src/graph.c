#include "graph.h"

/* create a new graph */
graph* new_graph()
{
    graph* g = ckallocz(sizeof(graph));
    g->node_hash = new_hashtable(8);
    return g;
}

/* create a new node for the graph */
node* new_node(const char* const name, void* const val)
{
    node* n = ckallocz(sizeof(node));
    
    n->name = copy_string(name);
    n->val  = val;

    n->key = -1;    
    n->visited = FALSE;    

    /* no edges initially */
    n->num_edges = 0;
    n->edges = NULL;

    return n;    
}

/* return TRUE if an edge exists from n1 to n2 */
static bool edge_exists(const node* const n1, const node* const n2)
{
    pre(n1 != NULL);
    pre(n2 != NULL);
    pre(n1 != n2);

    int i;
    for(i = 0; i < n1->num_edges; i++){
        forceassert(n1->edges[i]->n1 == n1 || n1->edges[i]->n2 == n1);
        if(n1->edges[i]->n1 == n2 || n1->edges[i]->n2 == n2){
            return TRUE;
        }
    }

    return FALSE;
}

/* make an edge from n1 to n2 */
static void make_edge(graph* const g, 
                      node* const n1, 
                      node* const n2, 
                      const int weight)
{
    /* no self edges */
    pre(n1 != n2);
    pre(n1 != NULL);
    pre(n2 != NULL);

    /* does an edge exist from n1 to n2 */
    if(!edge_exists(n1, n2)){
        /* create the edge. */
        edge* e    = ckalloc(sizeof(edge));
        e->next    = NULL;

        e->n1   = n1;
        e->n2   = n2;
        e->weight = weight;
        
        /* add it both of the nodes edge list */
        n1->num_edges++;
        n1->edges = ckrealloc(n1->edges, n1->num_edges * sizeof(edge*));
        n1->edges[n1->num_edges - 1] = e;

        n2->num_edges++;
        n2->edges = ckrealloc(n2->edges, n2->num_edges * sizeof(edge*));
        n2->edges[n2->num_edges - 1] = e;

        sladdhead(&g->edge_list, e);
    }
}

static inline int get_alignment_score(const node* const n1, 
                                      const node* const n2,
                                      int* const pmax1,
                                      int* const pmax2,
                                      int* const pmm)
{
    pre(n1 != NULL);
    pre(n2 != NULL);

    return get_read_alignment_score(n1->val, n2->val, pmax1, pmax2, pmm);
}

// static inline void remove_edges(graph* const g, node* const n1)
// {
//     pre(n1 != NULL);
// 
//     int i;
//     node* n2;
// 
//     for(i = 0; i < n1->num_edges; i++){
//         n2 = n1->edges[i]->n1 == n1 ? n1->edges[i]->n2 : n1->edges[i]->n1;
//         forceassert(n2->edges[n2->num_edges - 1] == n1->edges[i]);
//         slremove(&g->edge_list, n1->edges[i]);
//         ckfree(n1->edges[i]);
//         n2->num_edges--;
//         if(0 == n2->num_edges){
//             n2->edges = NULL;
//         }else{
//             n2->edges = ckrealloc(n2->edges, n2->num_edges * sizeof(edge*));
//         }
//     }
// 
//     n1->num_edges = 0;   
//     return;
// }

#if 0
/* add this node to the graph. Return the first node on the list that aligns
 * properly with this node. That will act as the input for the next node to be
 * added. This is only possible because the input is sorted */
node* add_node(graph* g, node* const n, node* const list)
{
    pre(g != NULL);
    pre(n != NULL);

    /* add this node to the hashtable of nodes. */
    bin* bin = add_hashtable(g->node_hash, n->name, strlen(n->name), n);
    ckfree(bin->name);
    *(&bin->name) = n->name;
    
    seqread* r1 = n->val;
    
    /* align it to other prospective nodes and create an edge if it passes
     * the thresholds */
    int score;
    seqread* r2;
    node* iter = list;

    /* the first node this node aligns to */
    node* probablecandidate = NULL;

    /* roughly the number of bases that align from read r1 and r2, and the
     * number of mismatches in the alignement */
    int max1, max2, diffs;

    while(iter){
        r2 = iter->val;
        assert(r1->s >= r2->s);
    
        if(r1->s >= r2->s && r1->s < r2->e){
            score = get_alignment_score(n, iter, &max1, &max2, &diffs);

            /* for considering two reads to be aligned, the following conditions
             * need to be satisfied: 
             * a) Overlap score needs to exceeds a threshold 
             * b) The overlap length should be >= than their speculated overlap
             *    on the reference.
             * c) Their speculated overlap on the reference should be should be
             *    greater than 5 bp.
             * d) The number of differences should be less than 10% of the 
             *    length of the overlap .*/
            if(score > THRESHOLD                    && 
               MIN(max1,max2) >= abs(r2->e - r1->s) &&
               abs(r2->e - r1->s) > 5               && 
               diffs <= (0.1 * MAX(max1, max2))){
                make_edge(g, iter, n, score);
                if(probablecandidate == NULL) probablecandidate = iter;
            }
        }

        if(r1->s > r2->e){
            break;
        }
        iter = iter->next;
    }

    /* now add this to the node list */ 
    sladdhead(&g->node_list, n);

    if(probablecandidate == NULL) probablecandidate = g->node_list;
    return probablecandidate;
}
#endif

/* add this node to the graph. Return the first node on the list that aligns
 * properly with this node. That will act as the input for the next node to be
 * added. This is only possible because the input is sorted */
void add_node(graph* g, node* const n)
{
    pre(g != NULL);
    pre(n != NULL);

    /* add this node to the hashtable of nodes. */
    bin* bin = add_hashtable(g->node_hash, n->name, strlen(n->name), n);
    ckfree(bin->name);
    *(&bin->name) = n->name;
    
    seqread* r1 = n->val;
    seqread* r2 = NULL;
    
    /* align it to other prospective nodes and create an edge if it passes
     * the thresholds */
    node* iter = g->node_list;

    /* roughly the number of bases that align from read r1 and r2, and the
     * number of mismatches in the alignement */
    int score, max1, max2, diffs, expectedoverlap;

    /* make edges with these nodes */
    node** nodes = NULL;
    int*  scores = NULL;
    int numnodes = 0;

    /* if a read aligns with one read, then it should align with all the reads
     * from then on. This should happen because the reads are sorted. If this 
     * does not happen then maybe this read should not be included in the
     * assembly  */
    while(iter){
        r2 = iter->val;
        forceassert(r1->s >= r2->s);
    
        if(r1->s >= r2->s && r1->s < r2->e){
            /* the expected overlap between these two reads */
            expectedoverlap = r2->e - r1->s;

            /* the observed alignment */
            score = get_alignment_score(n, iter, &max1, &max2, &diffs);

            /* if the observed alignment/overlap close to the expected
             * alignment/overlap */
            if(expectedoverlap > 5 && 
               score > THRESHOLD   && 
               abs(MIN(max1, max2) - expectedoverlap) < 5){
                numnodes++;

                scores = ckrealloc(scores, numnodes * sizeof(int));
                nodes  = ckrealloc(nodes,  numnodes * sizeof(node*));

                scores[numnodes-1] = score;
                nodes [numnodes-1] = iter;
            }
        }

        iter = iter->next;
    }

    /* make the selected edges */   
    int i;
    for(i = 0; i < numnodes; i++){
        make_edge(g, nodes[i], n, scores[i]);
    }

    ckfree(nodes);
    ckfree(scores);

    /* now add this to the node list */ 
    sladdhead(&g->node_list, n);
}

static void print_bin(bin* b)
{
    pre(b != NULL);

    printf("\"%s\";\n", b->name);
    node* n = b->val;

    int i;
    for(i = 0; i < n->num_edges; i++){
        printf("\"%s\" -- \"%s\"\n", 
              n->edges[i]->n1->name, n->edges[i]->n2->name);
    }
}

/* print a graph in the GRAPHVIZ format */
void print_graph(const graph* const g)
{
    pre(g != NULL);

    printf("graph g {\n");
    printf("node [shape = circle];\n");
    func_hashtable(g->node_hash, print_bin);
    printf("}\n");
}

/* a function which returns TRUE if all the nodes of the graph are visited */
bool allnodes_visited(const graph* const g)
{
    pre(g != NULL);

    node* iter;
    for(iter = g->node_list; iter; iter = iter->next){
        if(iter->visited == FALSE) return FALSE;
    }
    return TRUE;
}

/* does a node by this identifier exist in the graph */
node* node_exists(const graph* const g, const char* const name)
{
    pre(g != NULL);
    pre(name != NULL);

    return find_value_hashtable(g->node_hash, name, strlen(name));
}

/* compare function  used to sort the items for Prims Algorithm*/
static int compare(const void* const a, const void* const b)
{
    node* n1 = *((node**)a);
    node* n2 = *((node**)b);

    if(n2->key == n1->key){
        seqread* r1 = n1->val;
        seqread* r2 = n2->val;

        return r1->s - r2->s;
    }

    return n2->key - n1->key;
}

/* sort the edges based on the weights */
static int edgecomp(const void* const a, const void* const b)
{
    return (*((edge**)b))->weight - (*((edge**)a))->weight;
}

/* extract the node with the highest key */
static node* extract_max(node** nodes)
{
    node* n  = NULL;

    if((*nodes)->next != NULL){
        n = slpop(nodes);   
    }else{
        n = *nodes;
        *nodes = NULL;
    }
    forceassert(n != NULL);

    n->visited = TRUE;

    return n;
}

static void updateQ(node** const pQ, const node* const u, const int updated)
{
    /* iterate through Q till you find 'updated' elements that are connected to
     * 'u'. */
    int numfound = 0;

    node* iter = *pQ;
    for(; iter; iter = iter->next){
        if(edge_exists(u, iter)) numfound++;
        if(numfound == updated) break;
    }
    
    if(iter != NULL && iter->next != NULL && iter->next->next != NULL){
        node* tmp1 = iter->next;
        node* tmp2 = tmp1->next;

        tmp1->next = NULL;
        slsort(pQ, compare);

        if(tmp1->next == NULL) tmp1->next = tmp2;
        else{
            tmp1 = sllast(*pQ);
            tmp1->next = tmp2;
        }
    }else{
        slsort(pQ, compare);
    }
}

#ifndef NDEBUG
static bool issorted(const node* const list)
{
    const node* iter = list;
    for(; iter && iter->next; iter = iter->next){
        if(iter->key < iter->next->key){
            return FALSE;
        }
    }
    return TRUE;
}
#endif

/* find the maximum spanning tree for the given graph. Return a list of edges
 * which form a part of that path. The variables used in this routine are from
 * the Prim algorithm in "Introduction to Algorithms", by Thomas, Charles,
 * Ronald.  */
edge** find_maximum_spanning_tree(graph* const g, int* const numedges)
{
    pre(g != NULL);
    
    /* the path is empty initially*/
    *numedges = 0;
    edge** path = NULL; 
    node* list  = NULL;
    node* Q = g->node_list;
    
    /* set the key of the first node to be the maximum */
    Q->key = 0;

    int i, updated;
    node* t;
    node* u;
    node* v = NULL;
    while(Q){
        /* extract the next node */
        u = extract_max(&Q);
        sladdhead(&list, u);

        /* is there a path we can add to the spanning tree ? */
        qsort(u->edges, u->num_edges, sizeof(edge*), edgecomp);
        for(i = 0; i < u->num_edges; i++){
            t = u->edges[i]->n1 == u ? u->edges[i]->n2 : u->edges[i]->n1;
            if(t->visited == TRUE){ 
                *numedges += 1;
                path = ckrealloc(path, (*numedges) * sizeof(edge*));
                path[(*numedges) - 1] = u->edges[i];
                break;
            }
        }

        updated = 0;
        for(i = 0; i < u->num_edges; i++){  
            v = u->edges[i]->n1 == u ? u->edges[i]->n2 : u->edges[i]->n1;

            if(v->visited == FALSE){
                updated++;
                if(u->edges[i]->weight > v->key){
                    v->key = u->edges[i]->weight;
                }
            }
        }

        if(u->num_edges > 0) updateQ(&Q, u, updated);
        assert(issorted(Q));
        v = u;
    }

    g->node_list = list;

    return path;
}


/* free the resources held by this graph */
void free_graph(graph** pg)
{
    graph* g = *pg;
    
    free_hashtable(&g->node_hash);

    node* iter = g->node_list;
    for(; iter; iter = iter->next){
        ckfree(iter->edges);
    }

    slfreelist(&g->node_list);  
    slfreelist(&g->edge_list);

    ckfree(g);
}

/* count and return the number of singleton nodes in this graph */
uint countsingletons(const graph* const g)
{
    pre(g != NULL);

    uint count = 0;
    node* iter = g->node_list;
    for(; iter; iter = iter->next){
        if(iter->num_edges == 0) count++;    
    }

    return count;
}
