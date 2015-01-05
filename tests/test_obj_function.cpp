#include "my_utils.h"

using namespace std;

// ricordarsi che gli indici di riga e colonna in c++ partono da 0 !!

int main(int argc,char **argv)
{

	// mat A = mat("1 1 1 1; 0 0 1 1; 1 1 1 1");
	// cout << norm (A,"fro") << endl;
	// cout << norm (A,"fro")*norm (A,"fro")<< endl;
	// cout << pow(norm (A,"fro"),2) << endl;
	// cout << A << endl;
	// cout << "Numero righe:" << A.n_rows << endl;
	// cout << "Numero colonne:" << A.n_cols << endl;
	// cout << abs(-1) << endl;
	// size_t t(5);
	// cout << t << endl;

	mat U(4,4,fill::eye);
	mat H(4,4,fill::eye);
	sp_mat V(4,4);
	V(0,0) = V(1,1) = V(2,2) = V(3,3) = 1;

	sp_mat URM(4,4);
	URM(0,0) = 3;

	double lambda = 0;

	(U*H*V).print("UHV");

	URM.print("URM");

	double res =  evaluate_Obj_Function(URM,U,H,V,U,H,V,lambda);

	cout <<( (res ==4)? "Success" : "Failed" ) << endl;

	lambda = 1;

	res =  evaluate_Obj_Function(URM,U,H,V,U,H,V,lambda);

	cout <<( (res ==12)? "Success" : "Failed" ) << endl;

	return 0;
}
