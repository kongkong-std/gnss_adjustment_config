#include "../include/graph_impl.h"

/*
 * add node to adjacency list
 * */
AdjListNode * AdjListAddNode( int dest, double * weight, int weight_size )
{
    AdjListNode * new_node = NULL;
    if( ( new_node = ( AdjListNode * ) malloc( sizeof( AdjListNode ) ) )
            == NULL )
    {
        fprintf( stderr, "Memory allocation failed!\n" );
        exit( EXIT_FAILURE );
    }

    new_node->dest = dest;
    for( int index = 0; index < weight_size; index++ )
    {
        new_node->weight[ index ] = *( weight + index );
    }
    new_node->next = NULL;

    return new_node;
}

/*
 * graph generation with n vertices
 * */
AdjGraph * GraphGeneration( int n )
{
    AdjGraph * graph = NULL;
    if( ( graph = ( AdjGraph * ) malloc( sizeof( AdjGraph ) ) )
            == NULL )
    {
        fprintf( stderr, "Memory allocation failed!\n" );
        exit( EXIT_FAILURE );
    }

    graph->count_vertex = n;
    graph->count_edge = 0;

    if( ( graph->array = ( AdjList * ) malloc( graph->count_vertex * sizeof( AdjList ) ) )
            == NULL )
    {
        fprintf( stderr, "Memory allocation failed!\n" );
        exit( EXIT_FAILURE );
    }

#if 0
    printf( "address of graph = %p\n", graph );
#endif

    return graph;
}

/*
 * add edge to graph
 *     location[0] -> location[1]
 * */
void GraphInsert( AdjGraph * graph, int * location, double * weight, int weight_size )
{
    /*
     * location[0] src vertex
     * location[1] dest vertex
     * edge: src -> dest
     * */
    AdjListNode * current_node = NULL;    // current node in adjacency list
    AdjListNode * new_node = AdjListAddNode( *( location + 1 ), weight, weight_size );

    // add edge src -> dest
    if( graph->array[ *location ].head == NULL )
    {
        // current adjacency list head is null
        /*
         * add new node to adjacency list
         * */
        new_node->next = graph->array[ *location ].head;
        graph->array[ *location ].head = new_node;
    }
    else
    {
        // add new node from tail of linked list
        current_node = graph->array[ *location ].head;
        while( current_node->next != NULL )
        {
            current_node = current_node->next;
        }
        current_node->next = new_node;
    }
}

/*
 * graph display
 * */
void GraphDisplay( AdjGraph * graph )
{
    printf( "count of vertices(vertex): %d\n", graph->count_vertex );
    printf( "count of edges(edge): %d\n", graph->count_edge );

    // print adjacency list
    for( int index = 0; index < graph->count_vertex; index++ )
    {
        AdjListNode * head_node = graph->array[ index ].head;
        printf( "\nadjacancy list of vertex %d\n %d ", index, index );
        while( head_node != NULL )
        {
            printf( "-> %d ", head_node->dest );
            head_node = head_node->next;
        }
        putchar( '\n' );
    }
}


/*
 * destroy graph, free memory
 * */
void GraphDestroy( AdjGraph * graph )
{
    for( int index = 0; index < graph->count_vertex; index++ )
    {
        AdjListNode * current = NULL;
        while( graph->array[ index ].head != NULL )
        {
            current = ( graph->array[ index ].head )->next;
            free( graph->array[ index ].head );
            graph->array[ index ].head = current;
        }
    }

    free( graph->array );
    free( graph );
}

