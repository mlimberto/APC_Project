//
//  amf.h
//  
//


#ifndef ____amf__
#define ____amf__

#include <iostream>
#include <armadillo>

#include "pg_u.h" 
#include "pg_h.h"

using namespace arma;


// Definition of class AMF (Augmented matrix factorization)

// This class contains the implementation of the
// Augmented matrix factorization algorithm
// described in the paper

class AMF
{
    // ATTRIBUTES
    
        // Parameters
        double lambda_ ; // Over-fitting parameter
        int n_max_iter_ ; // maximum number of iterations
        double toll_ ; // Stop criterium 
        int r_; // number of latent factors
        int n_; // users
        int m_; // items
        int k_; // labels
    
        // Matrices
        SpMat<double> URM_;  // user-rating matrix, given in the dataset
        SpMat<int> ICM_; // item-content matrix, given in the dataset

        SpMat<double> V_,V_old_ ; // Matrix V
        Mat<double> U_,U_old_,H_,H_old_ ; // Matrices U,H



    
public :
    
    AMF(); // Constructor
    
    
    
    void solve();
    
};

#endif /* defined(____amf__) */
