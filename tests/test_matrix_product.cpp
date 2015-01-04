/* * test_matrix_product.cpp
 *
 *  Created on: Jan 4, 2015
 *      Author: andrea
 */


#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


int main(int argc,char **argv)
{
	mat A = mat("1 1 1 1; 0 0 1 1; 1 1 1 1");
	rowvec B = rowvec("3 3 3");
	colvec C = colvec("2 3 2 3");
	cout << A << endl;
	cout << B << endl;
	cout << C << endl;

	// il risultato deve essere 75 ed è effettivamente così =)
	// vediamolo

	mat res = B*A*C;
	cout << res << endl;

	return 0;
}
