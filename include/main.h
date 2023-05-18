#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
typedef struct
{
    int cnt_Station;
    int cnt_BaseStation;
    int cnt_RoverStation;
    int id_BaseStation[MAXSTATION];
    int id_RoverStation[MAXSTATION];
    double coo_BaseStation[TRI_MAXSTATION];
    int cnt_BaseLine;
    int mtd_Weight;
} Config;

/*
 * solution struct
 * rover station id
 * rover station position
 * least square result norm || Ax - b ||
 */
typedef struct
{
    int id_RoverStation[MAXSTATION];
    double coo_RoverStation[TRI_MAXSTATION];
    double norm_result;
} Solution;

// function prototype
/*
 * process of configure file
 * get value of configure file to Config type variables
 * input: FILE *, char *, Config *
 * output: Config *
 */
void ProcessConfigureFile(FILE *, char *, Config *);
void ProcessConfigureFile_sort(FILE *, char *, Config *);

void free_conf(Config *); // free memory

/*
 * process of source data file
 * get value of source data file to AdjGraph type variable
 * input: FILE *, char *, Config *, AdjGraph *
 * output: AdjGraph *
 */
void ProcessSourceData(FILE *, char *, Config *, AdjGraph *);

/*
 * encode station id, start from 0
 * base station in priority, then rover station
 * input: Config *, int, int *
 * output: int *
 */
void EncodeStationId(Config *, int, int *);

/*
 * implementation of adjustment of gnss network
 * input: Config *, AdjGraph *
 */
void ImplNetworkAdjustment(Config *, AdjGraph *, Solution *);

#endif
