/**
 * @file graph_impl.c
 * @author Zikang Qin
 * @brief graph data structure with adjacency list
 * @version 0.1
 * @date 2023-06-21
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "../include/graph_impl.h"

void InitializeLinkedList(AdjList *);

/**
 * @callgraph
 * @brief information of added edge,
 * destinated node of edge and weight data
 * of edge, assigning values to the struct
 *
 * @param [in] dest destinated node of edge
 * @param [in] weight weight of edge
 * @param [in] weight_size size of weight data in edge
 * @return AdjListNode* added edge
 */
AdjListNode *AdjListAddNode(int dest, double *weight, int weight_size)
{
    AdjListNode *new_node = NULL;
    if ((new_node = (AdjListNode *)malloc(sizeof(AdjListNode))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    new_node->dest = dest;
    for (int index = 0; index < weight_size; index++)
    {
        new_node->weight[index] = *(weight + index);
    }
    new_node->next = NULL;

    return new_node;
}

/**
 * @callgraph
 * @brief graph generation with n vertices
 * and n linked lists
 *
 * @param [in] n vertices of a graph
 * @return AdjGraph* graph data structure
 */
AdjGraph *GraphGeneration(int n)
{
    AdjGraph *graph = NULL;
    if ((graph = (AdjGraph *)malloc(sizeof(AdjGraph))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    graph->count_vertex = n;
    graph->count_edge = 0;

    if ((graph->array = (AdjList *)malloc(graph->count_vertex * sizeof(AdjList))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

#if 0
    printf( "address of graph = %p\n", graph );
#endif

#if 1 // initialize linkedlist in graph
    for (int index = 0; index < graph->count_vertex; index++)
    {
        InitializeLinkedList(graph->array + index);
    }
#endif

    return graph;
}

/**
 * @callgraph
 * @brief initialize linked list,
 * the head pointer of linked list to NULL
 * 
 * @param [in,out] pList linked list
 */
void InitializeLinkedList(AdjList *pList)
{
    pList->head = NULL;
}

/**
 * @callgraph
 * @brief add edge to graph, the numbering of 
 * the two vertices of an edge and weight data 
 * of the edge. source node is the first element of
 * location array, destinated node is the second
 * element of the edge. check if the head pointer
 * of the lined list is NULL, if NULL, add edge to
 * new adjacency list, else, add edge to tail of
 * current adjacency list
 * 
 * @param [in,out] graph graph data structure
 * @param [in] location source vertex and destinated vertes
 * @param [in] weight weight data of edge
 * @param [in] weight_size size of weight data
 */
void GraphInsert(AdjGraph *graph, int *location, double *weight, int weight_size)
{
    /*
     * location[0] src vertex
     * location[1] dest vertex
     * edge: src -> dest
     * */
    AdjListNode *current_node = NULL; // current node in adjacency list
    AdjListNode *new_node = AdjListAddNode(*(location + 1), weight, weight_size);

    // add edge src -> dest
    if (graph->array[*location].head == NULL)
    {
        // current adjacency list head is null
        /*
         * add new node to adjacency list
         * */
        new_node->next = graph->array[*location].head;
        graph->array[*location].head = new_node;
    }
    else
    {
        // add new node from tail of linked list
        current_node = graph->array[*location].head;
        while (current_node->next != NULL)
        {
            current_node = current_node->next;
        }
        current_node->next = new_node;
    }
}

/**
 * @callgraph
 * @brief display the graph with adjacency list,
 * graph traversal by vertex
 * 
 * @param [in] graph graph structure
 */
void GraphDisplay(AdjGraph *graph)
{
    printf("count of vertices(vertex): %d\n", graph->count_vertex);
    printf("count of edges(edge): %d\n", graph->count_edge);

    // print adjacency list
    for (int index = 0; index < graph->count_vertex; index++)
    {
        AdjListNode *head_node = graph->array[index].head;
        printf("\nadjacancy list of vertex %d\n %d ", index, index);
        while (head_node != NULL)
        {
            printf("-> %d ", head_node->dest);
            head_node = head_node->next;
        }
        putchar('\n');
    }
}

/**
 * @callgraph
 * @brief free memory of graph data structure
 * 
 * @param [in,out] graph grph data structure
 */
void GraphDestroy(AdjGraph *graph)
{
    for (int index = 0; index < graph->count_vertex; index++)
    {
        AdjListNode *current = NULL;
        while (graph->array[index].head != NULL)
        {
            current = (graph->array[index].head)->next;
            free(graph->array[index].head);
            graph->array[index].head = current;
        }
    }

    free(graph->array);
    free(graph);
}
