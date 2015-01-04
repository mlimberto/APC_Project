//
//  my_utils.h  
//
//  Public functions that all classes can use


#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <armadillo> 

using namespace arma;


///////////////////////////////////
////// INPUT-OUTPUT FUNCTIONS /////
///////////////////////////////////




///////////////////////////////////
//////// UTILITY FUNCTIONS ////////
///////////////////////////////////

// This function computes S(i,j) from URM,U,H,V using
// S(i,j) = URM(i,j) if URM(i,j) != 0 
// or 
// S(i,j) = (U*H*V)(i,j) if URM(i,j) == 0
double build_S(uword i, uword j,const SpMat<double>& URM,const Mat<double>& U,const Mat<double>& H,const SpMat<double>& V ) ;

void get_Positive_Matrix(Mat<double> &U) ;


///////////////////////////////////
////// DEBUGGING FUNCTIONS ////////
///////////////////////////////////


// Returns true if constraints on V are satisfied
// X is the ICM matrix
bool check_V_Constraint(const SpMat<double> &V,const SpMat<int>& X) ; 

bool check_Positive_Matrix(const Mat<double> &U) ;



#endif /* defined(UTILS_H) */

