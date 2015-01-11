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

// This function evaluates ||S - UHV ||^2  + lambda * (|| U ||^2 + || H ||^2 ) where
// ||*|| refers to the Frobenius norm 



double evaluate_Obj_Function(const SpMat<double>& URM,const Mat<double>& U,
							 const Mat<double>& H,const SpMat<double>& V,
							 const Mat<double>& U_old,const Mat<double>& H_old,
							 const SpMat<double>& V_old, const double lambda);

void get_Positive_Matrix(Mat<double> &U) ;


///////////////////////////////////
////// DEBUGGING FUNCTIONS ////////
///////////////////////////////////


// Returns true if constraints on V are satisfied
// X is the ICM matrix
bool check_V_Constraint(const SpMat<double> &V,const sp_umat &X) ;

bool check_Positive_Matrix(const Mat<double> &U) ;



#endif /* defined(UTILS_H) */

