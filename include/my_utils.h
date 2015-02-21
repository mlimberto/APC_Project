//
//  my_utils.h  
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto
//
//  Useful functions that are not class specific



#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <armadillo> 


///////////////////////////////////
////// INPUT-OUTPUT FUNCTIONS /////
///////////////////////////////////

template<typename T>
bool import_Sparse_Matrix(std::string mfilename,arma::SpMat<T> &MM,arma::umat &location_mat,arma::Col<T> &values);

template<typename T>
bool import_Sparse_Matrix(std::string mfilename,arma::SpMat<T> &MM);


///////////////////////////////////
//////// UTILITY FUNCTIONS ////////
///////////////////////////////////

// This function computes S(i,j) from URM,U,H,V using
// S(i,j) = URM(i,j) if URM(i,j) != 0 
// or 
// S(i,j) = (U*H*V)(i,j) if URM(i,j) == 0
double build_S(arma::uword i, arma::uword j,const arma::sp_mat& URM,const arma::mat& U,const arma::mat& H,const arma::sp_mat& V ) ;

arma::mat build_S_by_column(arma::uword j,const arma::sp_mat& URM,const arma::mat& U,const arma::mat& H,const arma::sp_mat& V );

arma::mat project_URM(const arma::sp_mat& URM,const arma::mat &S);



// This function evaluates ||S - UHV ||^2  + lambda * (|| U ||^2 + || H ||^2 ) where
// ||*|| refers to the Frobenius norm 
double evaluate_Obj_Function(const arma::sp_mat& URM,const arma::mat& U,
							 const arma::mat& H,const arma::sp_mat& V,
							 const arma::mat& U_old,const arma::mat& H_old,
							 const arma::sp_mat& V_old, const double lambda);

// This functions projects a matrix on the space of positive matrices
void get_Positive_Matrix(arma::mat &U) ;

arma::uvec get_Vector_Of_Indices(const arma::umat &L);

///////////////////////////////////
////// DEBUGGING FUNCTIONS ////////
///////////////////////////////////


// Returns true if constraints on V are satisfied
// X is the ICM matrix
bool check_V_Constraint(const arma::sp_mat &V,const arma::sp_umat &X) ;

bool check_Positive_Matrix(const arma::mat &U) ;



#endif /* defined(UTILS_H) */

