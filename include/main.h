#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph_impl.h"

// macro definition
#define MAXBUFFER 256
#define MAXSTATION 100
#define TRI_MAXSTATION 3 * MAXSTATION

// struct definition
/*
* basic adjustment configure
* count of stations in network
* count of base stations
* count of rover stations
* list of base stations
* list of rover stations
* coordinate of base stations
* count of baselines
* weight method to adjustment of gnss network
*     0: equal weight, 1: diagnoal covariance matrix, 2: full covariance matrix
*/
typedef struct {
    int cnt_Station;
    int cnt_BaseStation;
    int cnt_RoverStation;
    int id_BaseStation[MAXSTATION];
    int id_RoverStation[MAXSTATION];
    double coo_BaseStation[TRI_MAXSTATION];
    int cnt_BaseLine;
    int mtd_Weight;
}Config;

// function prototype
/*
* process of configure file
* get value of configure file to Config type variables
*/
void ProcessConfigureFile( FILE *, char *, Config * );
void ProcessConfigureFile_sort( FILE *, char *, Config * );

void free_conf( Config * );    // free memory

/*
* process of source data file
* get value of source data file to AdjGraph type variable
*/
void ProcessSourceData( FILE *, char *, Config *, AdjGraph * );

/*
* implementation of adjustment of gnss network
*/
void ImplNetworkAdjustment( Config *, AdjGraph * );

#endif
