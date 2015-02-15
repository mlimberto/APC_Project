#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"
#include "amf.h"

using namespace arma;
using namespace std;

bool test_PG_H(AMF & amf )
{
    // Set H,V,U_old

    amf.H_old_.eye(3,3);
    amf.H_old_(2,1) = 1.0;
    amf.H_old_(0,2) = 1.0;

    sp_mat V(3,5);
    V(0,0) = V(1,1) = V(2,2) = 1.0;
    V(2,3) = 0.5;
    V(0,3) = 0.5;
    V(1,4) = 1.0;

    swap(V,amf.V_old_);

    amf.U_old_.eye(5,3);
    amf.U_old_ = amf.U_old_*10.0;

    amf.U_.print("U_ matrix");


    // Set URM

    sp_mat URM(5,5);
    URM(0,0) = URM(1,1) = URM(2,2) = URM(3,3) = URM(4,4) =3.0;
    URM(1,1) = 4.0;
    URM(3,2) = 2.0;

    URM.print("URM Matrix");

    amf.URM_Tr_ = URM;

    // Initialize U_

    amf.U_.randu(5,3);
    amf.U_ = amf.U_*10.0;

    amf.U_.print("U_ matrix");

    amf.V_old_.print("V matrix");


    // Run the test

    amf.solve_pg_H_With_Log();

    amf.H_.print("H matrix");

    cout << "Lambda = " << amf.lambda_ << endl;

    cout << "End of the test function" << endl;

    return true;
}

int main(int argc,char** argv)
{
    // Create an instance of the AMF solver

    std::cout << "Creating instance of AMF ... "  << std::endl;
    AMF amf;
    //AMF amf("../APC_Project/dataset/Sampled_Dataset/urm_sampling.txt",
    //    "../APC_Project/dataset/Sampled_Dataset/icm_sampling.txt",
    //    "../APC_Project/tests/test_parameters.txt");

    //amf.initialize_Parameters("test_parameters.txt") ;
    amf.set_lambda(2.5);
    amf.set_n_max_iter(10000);
    amf.set_toll(1e-3);
    amf.set_n_latent_factors(50);
    amf.set_n_max_iter_gradient(800);
    amf.set_toll_gradient(1e-5);
    amf.set_gradient_step(1e-3);

    test_PG_H(amf);
    cout << ((test_PG_H(amf))? "Success " : "Failed" ) << endl;


    return 0;
}
