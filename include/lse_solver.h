#ifndef LSE_SOLVER_H_
#define LSE_SOLVER_H_

#include "adjustment_impl.h"

// function prototype
/*
 * equal weight
 * assemble least square equation
 * input: int **, double **, int, int, int,
 * int *, double **, int, int,
 * double **, double *, int, int
 * output: double **, double *
 */
void LSEAssemble_0(int **, double **, int, int, int,
                   int *, double **, int, int,
                   double **, double *, int, int);

/*
 * diagonal weight
 * assemble least square equation
 * input: int **, double **, int, int, int,
 * int *, double **, int, int,
 * double **, double *, int, int
 * output: double **, double *
 */
void LSEAssemble_1(int **, double **, int, int, int,
                   int *, double **, int, int,
                   double **, double *, int, int);

/*
 * full weight
 * assemble least square equation
 * input: int **, double **, int, int, int,
 * int *, double **, int, int,
 * double **, double *, int, int
 * output: double **, double *
 */
void LSEAssemble_2(int **, double **, int, int, int,
                   int *, double **, int, int,
                   double **, double *, int, int);

/*
 * determine code whether base station
 * input: int, int *, int
 * output: 1 = is base station, 0 = is not base station
 */
int IsBaseStation(int, int *, int);

/*
 * givens rotation
 * input: double * *, double *, int, int
 * ouput: double * *, double *
 */
void LSEGivens(double **, double *, int, int);

/*
 * implementation of givens rotation
 */
void givens_rotation_impl(double, double, int, int, int, int,
                          double **, double *);

/*
 * gaussian elimination solves uptriangular linear system
 * input: double * *, double *, int, double *
 * output: double *
 */
void LSEGauss(double **, double *, int, double *);

/*
 * assemble covariance matrix
 * input: double * *, int, double * *
 * output: double * *
 */
void assemble_mat_Cov(double **, int, double **);

/*
 * updating linear system with covariance matrix
 * input: double * *, int, int, double * *, double *
 * output: double * *, double *
 */
void update_linsys_full_weight(double **, int, int, double **, double *);

#endif