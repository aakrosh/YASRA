#include "graph.h"

/*return a new empty graph*/
GRAPH* new_graph()
{	
	GRAPH* g = ckallocz(sizeof(GRAPH));
	g->node_list = NULL;
	g->node_hash = new_hash(2);
	g->edge_list = NULL;
	return g;
}

/*free the resources help by the graph*/
void free_graph(GRAPH** pgraph)
{
	GRAPH* graph = *pgraph;
	NODE* iter = NULL;
	for(iter = graph->node_list; iter; iter = iter->next){
		ckfree(iter->edges);
	}
	slfreelist(&graph->node_list);
	slfreelist(&graph->edge_list);
	free_hash(&graph->node_hash);
	ckfree(graph);
}

/*add a new node to the graph*/
NODE* add_node(GRAPH* const g, void* const val, char* const name)
{
	NODE* node = NULL;

	/*is the node already there?*/
	if((node = hash_find_val(g->node_hash, name)) != NULL){
		return node;
	}

	node = ckallocz(sizeof(NODE));
	node->name = name;
	node->visited = FALSE;
	node->val = val;
	node->num_edges = 0;
	node->num_alloced = 0;
	node->edges = NULL;
	node->key = -1;

	sladd_head(&g->node_list, node);
	hash_add(g->node_hash, node->name, node);
	return node;
}

/*does an edge already exist?*/
static bool edge_exists(const GRAPH* const g, const NODE* const a, const NODE* const b)
{
	NODE* node = hash_must_find_val(g->node_hash, a->name);
	int i = 0;
	EDGE* edge = NULL;
	for(i = 0; i < node->num_edges; i++){
		edge = node->edges[i];
		if((edge->a == b) || (edge->b == b)){
			return TRUE;
		}
	}
    return FALSE;
}

/*add node to the edge*/
static void add_new_edge(NODE* const node, EDGE* const edge)
{
	/*do i need more memory ? */
	if(node->num_edges == node->num_alloced){
		int diff = node->num_alloced;
		int old_size = node->num_alloced * sizeof(EDGE);
		node->num_alloced +=  NUM_EDGES;
		int new_size = node->num_alloced * sizeof(EDGE);
		node->edges = ckrealloc(node->edges, new_size);
		memset(node->edges+diff, 0, (new_size-old_size));
	}

    node->edges[node->num_edges++] = edge;
}

/*add an edge in the graph from node a to node b*/
void add_edge(GRAPH* const g, NODE* const a, NODE* const b, const int weight)
{
	/*does the edge already exist? */
	if(edge_exists(g, a, b)){
		return;
	}

	/*if it doesnt exist, we go and make the to and fro connection*/
	EDGE* edge = ckallocz(sizeof(EDGE));
	edge->a = a;
	edge->b = b;
	edge->weight = weight;
	sladd_head(&g->edge_list, edge);

	add_new_edge(a, edge);
	add_new_edge(b, edge);
}

/* take the node, align it with every other node in the graph and if the 
 * overlap passes the thresholds then add an edge in the graph with the 
 * overlap score as the weight for the edge, return the number of new edges
 * formed.*/
int process_node(GRAPH* const g, NODE* const node, const bool suffix, int*** const pV, int** const pF, int*** pI)
{
	NODE* iter = NULL;
	NODE* list = g->node_list;
	int score = -1;
	int num_edges = 0;
	int bt = -1, et = -1;
	int br = -1, er = -1;
	sscanf(READ_HEAD((READ*)node->val),"%*s %d %d ", &bt, &et);
	int check = suffix ? OVERLAP_THRESHOLD : contained_threshold(et-bt);

	for(iter = list; iter; iter = iter->next){
		if(iter == node){
			continue;
		}

		sscanf(READ_HEAD((READ*)iter->val),"%*s %d %d ", &br, &er);
		if(er <= bt){
			break;
		}

		if((score = score_reads(iter->val, node->val, suffix, pV, pF, pI)) >= check){
#if DEBUG
			fprintf(stderr,"Edge %s : %s Score: %d\n", iter->name,node->name,score);
#endif
			add_edge(g, iter, node, score);
			num_edges++;
		}
	}

	return num_edges;
}

/* compare function  used to sort the items for Prims Algorithm*/
static int compare(const void* const a, const void* const b)
{
	return (*((NODE**)b))->key - (*((NODE**)a))->key;
}

/* find the node with the highest value of the key (it is the first element
 * in the list), remove it from the node list and node hash */
static NODE* extract_max(GRAPH* const g, NODE** Q)
{
	NODE* node = NULL;
	if(slcount(*Q) > 1){
		node = slpop(Q);
	}else{
		node = *Q;
		*Q = NULL;
	}
	if(node == NULL){
		fatal("Empty Q");
	}

	hash_remove_el(g->node_hash, node->name);

	return node;
}

static NODE* pop_earliest(NODE** plist)
{
	NODE* list = *plist;
	READ* seq = NULL;
	NODE* iter = NULL;
	NODE* min = NULL;
	int min_index = INT_MAX;
	int bt = 0;

	for(iter = list; iter; iter = iter->next){
		seq = iter->val;
		sscanf(READ_HEAD(seq),"%*s %d %*d ", &bt);
		if(bt < min_index){
			min_index = bt;
			min = iter;
		}
	}

	slremove(plist, min);
	return min;
}

/* find the maximum spanning tree for the graph. We progressively align the
 * nodes and return a list of contigs. The implementation is as per the
 * pseudocode in "Introduction to Algorithms". The implementation uses an array
 * right now. Changing it to Fibanacci heaps would lead to an improvement in
 * performance to O(E lg V).*/
CONTIG* MST(GRAPH* const g, NODE* const root, INTERVAL** const pIntervals, FILE* const rf, int*** const pV, int** const pF, int*** const pI)
{
	CONTIG* contig_list = NULL;
	CONTIG* contig = new_contig();
	NODE* list = NULL;
	NODE* u1 = NULL;
	NODE* u2 = NULL;
	NODE* Q = g->node_list;
	int i = 0, et = 0, bt = 0;

	/*lets set the key of the root to 0*/
	root->key = 0;
	slsort(&Q, compare);

	while(Q){
		/* extract the next node */
		u1 = extract_max(g, &Q);
		/*is it a new contig?*/
		if(u1->key == -1){
			sladd_head(&Q, u1);
			hash_add(g->node_hash, u1->name, u1);
			u1 = pop_earliest(&Q);
			sscanf(((READ*)u1->val)->header,"%*s %d %d ", &bt, &et);
			sladd_head(&contig_list, contig);
			contig = new_contig();
#if DEBUG
			fprintf(stderr, ">Contig: %s\n", READ_HEAD((READ*)u1->val));
#endif
		}
		/* add the returned node to a list, so that we can assign all the
		 * nodes in the graph to node_list*/
		sladd_head(&list, u1);
		/* align this new node to the existing profile */
		align_read(contig, u1->val, pIntervals, rf, pV, pF, pI);

		for(i = 0; i < u1->num_edges; i++){
			u2 = (u1->edges[i]->a == u1) ? u1->edges[i]->b : u1->edges[i]->a;
			if(hash_lookup(g->node_hash, u2->name) && (u1->edges[i]->weight > u2->key)){
				u2->key = u1->edges[i]->weight;
			}
		}
		slsort(&Q, compare);
	}

	g->node_list = list;
	sladd_head(&contig_list, contig);

	return contig_list;
}
