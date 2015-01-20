#ifndef GRAPH_H
#define GRAPH_H

#include "utilities.h"
#include "common.h"
#include "sequences.h"
#include "slinklist.h"
#include "hashtable.h"
#include "contig.h"

typedef struct node_st
{    
    struct node_st* next;
    char* name;     /* identifier for the node */
    void* val;      /* the actual value for the node */

    int key;        /* useful as key/weight */
    bool visited;   /* useful for testing membership, visited etc. */

    /* edges to and from this node */
    int num_edges;
    struct edge_st** edges;
}node;

/* On a 64 bit machine, this should be 8 + 8 + 8 + sizeof(int) */
typedef struct edge_st
{
    struct edge_st* next;
    node* n1;
    node* n2;
    int weight;
}edge;

typedef struct graph_st
{   
    node* node_list;        /* list of all the nodes in this graph */
    edge* edge_list;        /* list of all the edges in this graph */
    hashtable* node_hash;   /* hashtable of all the nodes in this graph */    
}graph;

/* create a new graph */
graph* new_graph();

/* create a new node for the graph. */
node* new_node(const char* const name, void* const val);

/* add this node to the graph */
void add_node(graph* g, node* const n);

/* print a graph in the GRAPHVIZ format */
void print_graph(const graph* const g);

/* a function which returns TRUE if all the nodes of the graph are visited */
bool allnodes_visited(const graph* const g);

/* does a node by this identifier exist in the graph */
node* node_exists(const graph* const g, const char* const name);

/* return the maximum spanning tree for this graph */
edge** find_maximum_spanning_tree(graph* const g, int* const numedges);

/* free the resources held by this graph */
void free_graph(graph** pg);

/* count and return the number of singleton nodes in this graph */
uint countsingletons(const graph* const g);

#endif
