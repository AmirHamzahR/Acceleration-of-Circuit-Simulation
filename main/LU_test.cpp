#include <iostream>
#include "armadillo"
using namespace arma;
using namespace std;

int main(){

   mat A = randu<mat>(5,5);
   vec b = randu<vec>(5);

    A.print("A =");
    b.print("b =");

    // solve Ax = b
    vec x = solve(A,b);

    // print x
    x.print("x =");

    // find LU decomp of A, if needed, P is the permutation matrix
    mat L, U, P;
    lu(L,U,P,A);

    // print l
    L.print("L = ");

    // print U
    U.print("U = ");
    
    //solve using PLU decomposition
    vec x1 = solve(trimatu(U), solve(trimatl(L), P*b));

    // print PLU decomposition solution
    x1.print("PLU solution:");

    //Check that A = LU
    //(A-(P.t()*L*U)).print("Test of LU decomposition:");

    return 0;
}