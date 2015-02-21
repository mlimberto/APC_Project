#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"
#include "amf.h"

using namespace arma;
using namespace std;

bool test_top_N(AMF & amf )
{
    mat URM(5,5,fill::eye);
    mat U(5,5,fill::ones);
    mat H(5,5,fill::eye);
    mat V(5,5,fill::eye);

    U(0,1) = 3;
    U(1,2) = 3;
    U(2,3) = 3;
    U(3,4) = 3;
    U(4,0) = 3;

    U(0,2) = 4;
    U(1,3) = 4;
    U(2,4) = 4;
    U(3,0) = 4;
    U(4,1) = 4;

    amf.URM_Tr_ = 5*URM ;
    amf.U_old_ = U;
    amf.H_old_ = H;
    amf.V_old_ = V;

    amf.n_ = URM.n_rows;
    amf.m_ = URM.n_cols;

    cout << amf.n_ << endl;
    cout << amf.m_ << endl;

    amf.URM_Tr_.print("URM matrix");

    (U*H*V).print("Recommendations");

    umat result = amf.get_TopN_Recommendation(2);

    result.print("Result");

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
    amf.set_n_max_iter(10000);
    amf.set_n_latent_factors(50);
    amf.set_n_max_iter_gradient(1);
    amf.set_toll_gradient(1e-5);

    test_top_N(amf);


    return 0;
}
