#ifndef ADJUSTMENT_IMPL_H_
#define ADJUSTMENT_IMPL_H_

#include "main.h"

// function prototype
/*
* equal weight implementation
*/
void ImplNetworkAdjustment_0( Config *, AdjGraph * );

/*
* diagonal weight implementation
*/
void ImplNetworkAdjustment_1( Config *, AdjGraph * );

/*
* full weight implementation
*/
void ImplNetworkAdjustment_2( Config *, AdjGraph * );

#endif