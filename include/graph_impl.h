#ifndef GRAPH_IMPL_H_
#define GRAPH_IMPL_H_

// header
#include <stdio.h>
#include <stdlib.h>

// graph implementation with adjacency list
/*
 * adjacency list node
 *     dest vertex location
 *     weight data of edge
 *     next node
 * */
typedef struct adj_list_node{
    int dest;    // dest vertex location
    double weight[9];    // weight array
    struct adj_list_node * next;    // next node
}AdjListNode;

/*
 * adjacency list
 *     head node
 * */
typedef struct adj_list{
    AdjListNode * head;    // head node in adjacency list
}AdjList;

/*
 * graph represents with adjacnecy list
 *     count of vertex(vertices)
 *     count of edges
 *     adjacency list array
 * */
typedef struct adj_graph{
    int count_vertex;
    int count_edge;
    AdjList * array;
}AdjGraph;

// function prototype
/*
 * add new node
 * */
AdjListNode * AdjListAddNode( int, double *, int );

/*
 * add edge to graph
 *     src vertex -> dest vertex
 *     weight data and its size
 * */
void GraphInsert( AdjGraph *, int *, double *, int );

/*
 * destroy graph
 *     free graph
 * */
void GraphDestroy( AdjGraph * );

/*
 * graph generation with vertices
 * */
AdjGraph * GraphGeneration( int );

/*
 * graph display
 *     count of vertex
 *     count of edge
 *     adjacency list
 * */
void GraphDisplay( AdjGraph * );

#endif
