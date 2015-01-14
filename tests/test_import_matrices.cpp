#include <iostream> 
#include <string>
#include <vector>

#include "my_utils.h"

using namespace std;
using namespace arma;

int main(int argc,char** argv)
{

    std::string matrixICM_filename("../dataset/Sampled_Dataset/icm_sampling.txt");

    std::string matrixURM_filename("../dataset/Sampled_Dataset/urm_sampling.txt");

    sp_umat M;


    import_Sparse_Matrix<uword>(matrixICM_filename,M);

    M.print();

    SpMat<double> urm;

    import_Sparse_Matrix<double>(matrixURM_filename,urm);

    urm.print();


    return 0;
}
