//
//  amf.h
//  
//

#ifndef ____amf__
#define ____amf__

#include <iostream>
#include <armadillo>

#include "my_utils.h"
#include "pg_u.h" 
#include "pg_h.h"
#include "compute_v.h"

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
        int n_max_iter_gradient_ ; // Same as above for the projected gradient method
        double toll_gradient_;
        int r_; // number of latent factors
        int n_; // users
        int m_; // items
        int k_; // labels
    
        // Matrices
        SpMat<double> URM_;  // user-rating matrix, given in the dataset
        SpMat<int> ICM_; // item-content matrix, given in the dataset

        SpMat<double> V_,V_old_ ; // Matrix V
        Mat<double> U_,U_old_,H_,H_old_ ; // Matrices U,H

    // METHODS
        void initialize_matrices();

        void solve_one_iteration();

    
public :

    // CONSTRUCTORS
    
    AMF();
    
    // SOLVING METHODS
        
    void solve();

    // SET METHODS

    inline void set_lambda(const double &lambda) { lambda_ = lambda; }
    inline void set_n_max_iter(int n) { n_max_iter_ = n; }
    inline void set_toll(const double &toll) { toll_ = toll; }
    inline void set_n_latent_factors(int r) { r_ = r; }
    inline void set_n_max_iter_gradient(int n) { n_max_iter_gradient_ = n; }
    inline void set_toll_gradient(const double &toll) { toll_gradient_ = toll; }
    
    // GET METHODS

    inline double get_lambda() { return lambda_; }
    inline int get_n_max_iter() { return n_max_iter_; }
    inline double get_toll() { return toll_; }
    inline int get_n_latent_factors() { return r_; }
    inline int get_n_max_iter_gradient() { return n_max_iter_gradient_; }
    inline double get_toll_gradient() { return toll_gradient_; }
};

#endif /* defined(____amf__) */
