#ifndef ADJUSTMENT_IMPL_H_
#define ADJUSTMENT_IMPL_H_

#include "main.h"

// function prototype
/*
 * equal weight implementation
 * input: Config *, AdjGraph *
 */
void ImplNetworkAdjustment_0(Config *, AdjGraph *, Solution *);

/*
 * diagonal weight implementation
 * input: Config *, AdjGraph *
 */
void ImplNetworkAdjustment_1(Config *, AdjGraph *, Solution *);

/*
 * full weight implementation
 * input: Config *, AdjGraph *
 */
void ImplNetworkAdjustment_2(Config *, AdjGraph *, Solution *);

/*
 * solver for least square equation
 * equal weight
 * input: int **, double **, int, int, int, double **, double *, int, int, double *, double *
 * output: double **, double *, double *, double *
 */
void LSESolver_0(Config *, int **, double **, int, int, int,
                 double **, double *, int, int, double *, double *);

/*
 * solver for least square equation
 * diagonal weight
 * input: int **, double **, int, int, int, double **, double *, int, int, double *, double *
 * output: double **, double *, double *, double *
 */
void LSESolver_1(Config *, int **, double **, int, int, int,
                 double **, double *, int, int, double *, double *);

/*
 * solver for least square equation
 * full weight
 * input: int **, double **, int, int, int, double **, double *, int, int, double *, double *
 * output: double **, double *, double *, double *
 */
void LSESolver_2(Config *, int **, double **, int, int, int,
                 double **, double *, int, int, double *, double *);

/*
 * initialize least square equation
 * mat = 0, sol = 0, b = 0, residual = 0
 */
void InitLSE(double *, double *, double *, double **, int, int);

/*
* computing L2 norm of vector a
* input: double *, int
* output: result
*/
double lse_norm_2( double *, int );

/*
* updating lse solution
* input: Config *, double *, double *, int, int, Solution *
* output: Solution *
*/
void LSEUpdateSolution( Config *, double *, double *, int, int, Solution * );

#endif