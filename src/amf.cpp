//
//  amf.cpp
//  
//
//
//

#include "amf.h"

//////////////////////////
////  CONSTRUCTORS  //////
//////////////////////////

// Most basic constructor
AMF::AMF() :
lambda_(0),n_max_iter_(100),toll_(0.001)
{
    
}



//////////////////////////
////  METHODS       //////
//////////////////////////


// TODO La funzione non verifica che tutti i parametri siano stati 
// inizializzati, sarebbe da fare ma non prioritario
bool AMF::initialize_Parameters(std::string filename)
{
	#ifndef NDEBUG
	std::cout << "Initializing parameters of AMF class from file " << filename << std::endl;
	#endif

	// Define the available options through a std::map object
	static std::map<std::string,unsigned int> options_map; 
	options_map["LAMBDA"] = 0;
	options_map["N_MAX_ITER"] = 1;
	options_map["TOLL"] = 2;
	options_map["N_LATENT_FACTORS"] = 3;
	options_map["N_MAX_ITER_GRADIENT"] = 4;
	options_map["TOLL_GRADIENT"] = 5;


	// Open and read the parameters file 
	std::ifstream param_file(filename);
	std::string option;

	while (param_file >> option)
	{
		#ifndef NDEBUG
		std::cout << "Setting option " << option << std::endl;
		#endif

		// Check if the option is valid
		std::map<std::string,unsigned int>::const_iterator it = options_map.find(option);

		if (it == options_map.end())
		{
			std::cerr << "ERROR : Check consistency of file " << filename << std::endl;
			return false; 
		}

		// Do the mapping and assign the parameter
		unsigned int mapped_option = options_map[option] ;

		switch(mapped_option)
		{
			case 0 : // Next input is value of lambda
			double lambda ; param_file >> lambda ; set_lambda(lambda);
			break;

			case 1:
			unsigned int n_max_iter ; param_file >> n_max_iter ; 
			set_n_max_iter(n_max_iter) ;
			break;

			case 2: 
			double toll ; param_file >> toll ; set_toll(toll);
			break;

			case 3:
			unsigned int n_lf ; param_file >> n_lf ; 
			set_n_latent_factors(n_lf);
			break;

			case 4:
			unsigned int n_max_iter_gradient ; param_file >>n_max_iter_gradient ; 
			set_n_max_iter_gradient(n_max_iter_gradient) ;
			break;

			case 5: 
			double toll_g ; param_file >> toll_g ; set_toll_gradient(toll_g);
			break;

			default : 
			std::cout << "Running default option" << endl;
			break;

		}
	}

	return true;
}

bool AMF::inizialize_ICM_Locations(std::string matrix_filename){
    #ifndef NDEBUG
    std::cout << "Importing ICM from file " << matrix_filename << std::endl;
    #endif

    arma::umat Location_Matrix;
    arma::uvec Values;
    arma::umat RCi(2,1);
    unsigned int val;
    std::ifstream matrix_file(matrix_filename);

    // Read the file and build the Location Matrix and the Values vector
    unsigned int i=0;
    while (matrix_file >> RCi(0,0) >> RCi(1,0) >> val){
        Location_Matrix.insert_cols(i,RCi);
        Values.resize(i+1);
        Values(i)=val;
        i++;
    }

    arma::sp_umat icm(Location_Matrix,Values);

    // Con lo swap la variabile urm prende il posto di URM_ 
    // mentre la vecchia URM_ viene distrutta dal distruttore
    // alla fine dell'esecuzione della funzione
    std::swap(ICM_,icm);

    return true;
}

bool AMF::inizialize_URM_Locations(std::string matrix_filename){
    #ifndef NDEBUG
    std::cout << "Importing URM from file " << matrix_filename << std::endl;
    #endif

    arma::umat Location_Matrix;
    arma::vec Values;
    arma::umat RCi(2,1);
    double val;
    std::ifstream matrix_file(matrix_filename);

    // Read the file and build the Location Matrix and the Values vector
    unsigned int i=0;
    while (matrix_file >> RCi(0,0) >> RCi(1,0) >> val){
        Location_Matrix.insert_cols(i,RCi);
        Values.resize(i+1);
        Values(i)=val;
        i++;
    }

    #ifndef NDEBUG
    Location_Matrix.print("Location Matrix : ");
    Values.print("Values Vector: ");
    #endif

    arma::sp_mat urm(Location_Matrix,Values);

    std::swap(URM_,urm);

    return true;
}

void AMF::initialize_matrices()
{

}


void AMF::solve() 
{

	// Remember to swap the variables

}

