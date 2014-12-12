//
//  my_utils.h  
//
//  Public functions that all classes can use


#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <armadillo> 

using namespace arma;


// Utility functions

// This function computes S(i,j) from URM,U,H,V using
// S(i,j) = URM(i,j) if URM(i,j) != 0 
// or 
// S(i,j) = (U*H*V)(i,j) if URM(i,j) == 0
double build_S(int i, int j, SpMat<double> URM, Mat<double> U, Mat<double> H, SpMat<double> V ) ; 


// Debugging functions

bool check_V_Constraint(SpMat<double> V) ; // Verify that constraints on V are satisfied



#endif /* defined(UTILS_H) */

