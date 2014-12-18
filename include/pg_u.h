//
//  pg_u.h
//  
//  Class for the projected gradient solver
//  for the matrix U

#ifndef PG_U_H
#define PG_U_H

#include <iostream>
#include <armadillo>

#include "my_utils.h"

using namespace arma;



class PG_U {




public :

	PG_U(); // Constructor

	void solve();


};









#endif /* defined(PG_U_H) */