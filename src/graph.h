#ifndef GRAPH_H
#define GRAPH_H

#include "hash.h"
#include "slist.h"
#include "contig.h"

#define NUM_EDGES 20

typedef struct graph
{
	struct node* node_list;
	struct edge* edge_list;
	struct hash* node_hash;
}GRAPH;

typedef struct node
{
	struct node* next;		/*the next in the list*/
	char* name;				/*the name of the node*/
	void* val;				/*the data in the node*/
	struct edge** edges;	/*edges from this node*/
	int num_edges;			/*number of edges*/
	int num_alloced;		/*number of edges allocated for*/
	int key;				/*the basis for sorting or for MST*/
	bool visited;			/*did we visit the node*/
}NODE;

typedef struct edge
{
	struct edge* next;		/*next edge in the master edge list*/
	struct node* a;			/*one of the nodes connected by the edge*/
	struct node* b;			/*the other node*/
	int weight;				/*weight of the edge*/
}EDGE;

/*return a new empty graph*/
GRAPH* new_graph();

/*free the resources help by the graph*/
void free_graph(GRAPH** pgraph);

/*add a new node to the graph*/
NODE* add_node(GRAPH* const g, void* const val, char* const name);

/* take the node, align it with every other node in the graph and if the 
 * overlap passes the thresholds then add an edge in the graph with the 
 * overlap score as the weight for the edge, return the number of new edges
 * formed.*/
int process_node(GRAPH* const g, NODE* const node, const bool suffix, int*** const pV, int** const pF, int*** const pI);

/*run Prims Algorithm, return the list of contigs which are the result*/
CONTIG* MST(GRAPH* g, NODE* const root, INTERVAL** const pIntervals, FILE* const rf, int*** const pV, int** const pF, int*** const pI);
#endif
