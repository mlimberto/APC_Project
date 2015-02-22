//
//  amf.h
//  
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto

#ifndef ____amf__
#define ____amf__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <armadillo>

#include "my_utils.h"

#ifdef AMFTIME
#include <chrono>
#endif

#ifdef AMFOPENMP
#include <omp.h>
#endif


// Definition of class AMF (Augmented matrix factorization)

// This class contains the implementation of the
// Augmented matrix factorization algorithm
// described in the paper

class AMF
{
    // ATTRIBUTES

        // Log file
        std::ofstream total_logfile_;

        std::string amf_filename_prefix_;
    
        // Parameters
        double lambda_ ; // Over-fitting parameter

        unsigned int n_max_iter_ ; // maximum number of iterations

        unsigned int n_max_iter_gradient_ ; // Same as above for the projected gradient method
        double toll_gradient_;
        double gradient_step_; // Length of the gradient step

        arma::uword r_; // number of latent factors
        arma::uword n_; // users
        arma::uword m_; // items
        arma::uword k_; // labels

    
        // Matrices and their coordinates
        arma::sp_mat URM_Tr_;  // user-rating matrix, given in the dataset
        arma::sp_mat URM_Val_;
        arma::sp_umat ICM_; // item-content matrix, given in the dataset


        arma::umat URM_Tr_Location_Matrix_;
        arma::umat URM_Val_Location_Matrix_;
        arma::umat ICM_Location_Matrix_;

        arma::vec URM_Val_Values_;
        arma::vec URM_Tr_Values_;
        arma::uvec ICM_Values_;

        arma::sp_mat V_,V_old_ ; // Matrix V
        arma::mat U_,U_old_,H_,H_old_ ; // Matrices U,H

    // METHODS
        void initialize_matrices();

        void solve_one_iteration();

        void solve_pg_U();

        void solve_pg_U_With_Log();

        void solve_pg_U_One_Iteration(arma::mat G,const arma::mat &A, const arma::mat &AAt);


        void solve_pg_H();

        void solve_pg_H_With_Log();

        void solve_pg_H_One_Iteration(arma::mat G,const arma::mat &UtU, 
                                        const arma::mat &VVt);


        void solve_V();

        void solve_V_With_Log();

        void solve_V_One_Iteration(arma::mat G,const arma::mat &WtW, const double prec_obj,double &curr_obj);

        void project_V(arma::sp_mat &V_new ,arma::mat &V_hat);    

        arma::sp_mat project_ICM(const arma::mat &G);

        double evaluate_Against_URM_Validation();
    
public :

    // CONSTRUCTORS and INITIALIZERS
    
    AMF();

    AMF(std::string URM_Tr_filename, std::string ICM_filename, std::string param_filename);

    bool initialize_Parameters(std::string filename);
    
    // SOLVING METHODS
        
    void solve();
    
    void solve_With_Log();

    void solve_For_Tuning(std::string mfilename);

    void solve_For_Tuning_With_Log(std::string mfilename);


    // RECOMMENDATIONS

    arma::umat get_TopN_Recommendation(arma::uword N, bool export_to_file = true);


    // SET METHODS

    inline void set_lambda(double lambda) { lambda_ = lambda; }
    inline void set_n_max_iter(unsigned int n) { n_max_iter_ = n; }
    inline void set_n_latent_factors(arma::uword r) { r_ = r; }
    inline void set_n_max_iter_gradient(unsigned int n) { n_max_iter_gradient_ = n; }
    inline void set_toll_gradient(double toll) { toll_gradient_ = toll; }
    
    // GET METHODS

    inline double get_lambda() { return lambda_; }
    inline unsigned int get_n_max_iter() { return n_max_iter_; }
    inline arma::uword get_n_latent_factors() { return r_; }
    inline unsigned int get_n_max_iter_gradient() { return n_max_iter_gradient_; }
    inline double get_toll_gradient() { return toll_gradient_; }
    inline double get_gradient_step() { return gradient_step_; }

    inline void print_ICM(){ICM_.print("ICM =");}
    inline void print_URM_Tr(){URM_Tr_.print("URM =");}

    // EXPORT TO FILE METHODS

    void export_Results();

    void import_Results();

};

#endif /* defined(____amf__) */
