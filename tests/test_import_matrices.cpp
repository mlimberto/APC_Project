#include <iostream> 
#include <string>
#include <vector>

#include "my_utils.h"
#include "amf.h"
int main(int argc,char** argv)
{


    // Create instance for class AMF

    AMF solver;


    std::string matrixICM_filename("../APC_Project/dataset/Matrices/icm_sampling.txt");
    solver.inizialize_ICM_Locations(matrixICM_filename);
    solver.print_ICM();
    std::string matrixURM_filename("../APC_Project/dataset/Matrices/urm_converted.txt");
    solver.inizialize_URM_Locations(matrixURM_filename);
    solver.print_URM();


    return 0;

    return 0;
}
