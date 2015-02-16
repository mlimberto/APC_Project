#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"
#include "amf.h"

using namespace arma;
using namespace std;

bool test_PG_V(AMF & amf )
{
    // Set URM
    sp_mat URM(3,3);
    URM(0,0) = URM(1,1) = URM(2,2) = 4.0;
    URM(0,1) = 1.0;
    URM(0,2) = 1.0;
    URM(1,0) = 1.0;
    URM(1,2) = 1.0;
    URM(2,0) = 1.0;
    URM(2,1) = 1.0;

    URM.print("URM Matrix");

    amf.URM_Tr_ = URM;

    // Set ICM
    umat locations ; 
    locations << 0 << 0 << 0 << 1 << 2 << endr
              << 0 << 1 << 2 << 1 << 2 << endr;
    uvec values ;
    values << 1 << 1 << 1 << 1 << 1 << endr;

    sp_umat ICM(locations,values);


    ICM.print("ICM matrix");


    amf.ICM_ = ICM;   

    amf.ICM_Location_Matrix_ = locations;
    amf.ICM_Values_ = values;

    // Initialize matrices
    amf.initialize_matrices(); 

    // Set H,V,U_
    amf.U_old_.eye(3,3);
    amf.U_old_ = 1.0*amf.U_old_;

    amf.H_old_.eye(3,3);
    amf.H_old_ = 1.0*amf.H_old_;

    amf.U_ = amf.U_old_ ;
    amf.H_ = amf.H_old_ ;

    amf.V_old_.print("V initialized");

    // sp_mat V(3,3);
    // V(0,0) = V(1,1) = V(2,2) = 1.0;

    // swap(V,amf.V_old_);

    // Run the test

    amf.solve_V_With_Log();

    amf.V_.print("V matrix");

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
    amf.set_lambda(1.0);
    amf.set_n_max_iter(10);
    amf.set_toll(1e-3);
    amf.set_n_latent_factors(50);
    amf.set_n_max_iter_gradient(1000);
    amf.set_toll_gradient(1e-5);
    amf.set_gradient_step(0.01);

    cout << ((test_PG_V(amf))? "Success " : "Failed" ) << endl;


    return 0;
}
