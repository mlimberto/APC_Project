/* * test_matrix_product.cpp
 *
 *  Created on: Jan 4, 2015
 *      Author: andrea
 */





#include "my_utils.h"

using namespace std;

// ricordarsi che gli indici di riga e colonna in c++ partono da 0 !!

int main(int argc,char **argv)
{

	mat A = mat("1 1 1 1; 0 0 1 1; 1 1 1 1");
	mat U = mat("3 3 3; 1 1 1; 1 1 1; 1 1 1");
	sp_mat V = sp_mat("2 0 0 0; 3 0 0 0; 2 0 0 0; 3 0 0 0");
	sp_mat URM = sp_mat("0 1 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0");


	cout << URM << endl;

	cout << U*A*V << endl;

	// NOTA: per le matrici sparse si accede e si cambiano i valori
	//       come per le matrici piene !! Fa tutto in automatico !
	// OTTIMA NOTIZIA

	/*URM(0,0)=5;
	cout << URM << endl;
	*/
	// cout << (U.row(i))*A*(V.col(j)) << endl;

	// il risultato deve essere 75 ed è effettivamente così =)
	// vediamolo

	cout <<"S(0,0) = " << build_S(0, 0, URM, U, A , V ) << endl;
	cout <<"S(0,1) = " << build_S(0, 1, URM, U, A , V ) << endl;

	return 0;
}
