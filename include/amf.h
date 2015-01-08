//
//  amf.h
//  
//

#ifndef ____amf__
#define ____amf__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <armadillo>

#include "my_utils.h"

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
        double gradient_step_; // Length of the gradient step

        int r_; // number of latent factors
        int n_; // users
        int m_; // items
        int k_; // labels
    
        // Matrices
        sp_mat URM_;  // user-rating matrix, given in the dataset
        sp_umat ICM_; // item-content matrix, given in the dataset

        SpMat<double> V_,V_old_ ; // Matrix V
        Mat<double> U_,U_old_,H_,H_old_ ; // Matrices U,H

    // METHODS
        void initialize_matrices();

        void solve_one_iteration();


        void solve_pg_U();

        void solve_pg_U_With_Log();

        void solve_pg_U_One_Iteration(mat &G,const mat &A);


        void solve_pg_H();

        void solve_pg_H_With_Log();

        void solve_pg_H_One_Iteration(mat &G, mat &A);


        void solve_V();

        sp_mat solve_V_One_Step_Gradient(const sp_mat &V_0);
        sp_mat project_ICM(const mat &G);

        void othogonal_projection(sp_mat &v);
    
public :

    friend bool test_PG_U(AMF & amf); // For testing purposes

    // CONSTRUCTORS and INITIALIZERS
    
    AMF();

    AMF(std::string URM_filename, std::string ICM_filename, std::string param_filename);

    bool initialize_Parameters(std::string filename);
    bool initialize_ICM_Locations(std::string matrix_filename);
    bool initialize_URM_Locations(std::string matrix_filename);
    
    // SOLVING METHODS
        
    void solve();

    // SET METHODS

    inline void set_lambda(double lambda) { lambda_ = lambda; }
    inline void set_n_max_iter(int n) { n_max_iter_ = n; }
    inline void set_toll(double toll) { toll_ = toll; }
    inline void set_n_latent_factors(int r) { r_ = r; }
    inline void set_n_max_iter_gradient(int n) { n_max_iter_gradient_ = n; }
    inline void set_toll_gradient(double toll) { toll_gradient_ = toll; }
    inline void set_gradient_step(double g_step) { gradient_step_ = g_step; }
    
    // GET METHODS

    inline double get_lambda() { return lambda_; }
    inline int get_n_max_iter() { return n_max_iter_; }
    inline double get_toll() { return toll_; }
    inline int get_n_latent_factors() { return r_; }
    inline int get_n_max_iter_gradient() { return n_max_iter_gradient_; }
    inline double get_toll_gradient() { return toll_gradient_; }
    inline double get_gradient_step() { return gradient_step_; }


    inline void print_ICM(){ICM_.print("ICM =");}
    inline void print_URM(){URM_.print("URM =");}

};

#endif /* defined(____amf__) */
