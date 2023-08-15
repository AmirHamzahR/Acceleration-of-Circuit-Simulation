#include <iostream>
#include <armadillo>
#include <tuple>
#include <algorithm>
#include <list>
using namespace arma;
using namespace std;

mat branch_ext(mat M, int node_x, int node_y){
    int M_cols = M.n_cols;
    int M_rows = M.n_rows;
    mat va = zeros(M_rows,1);
    if(node_x == 0){
        va.row(node_y-1).col(0) = 1;
    }
    else if(node_y == 0){
        va.row(node_x-1).col(0) = 1;
    }
    else{
        va.row(node_x-1).col(0) = -1;
        va.row(node_y-1).col(0) = 1;
    }

    mat zero_ext = zeros(1,1);
    mat ha = va.as_row();
    mat haz = join_rows(ha,zero_ext);

    mat M1 = join_rows(M,va);
    mat M2 = join_cols(M1,haz);

    return M2;
}

// mat mat_sum(mat list_of_matrices){
//     mat a = list_of_matrices[0];
//     for (int i = 0; i < list_of_matrices.n_elem; ++i) {
//         a = a + list_of_matrices[i];
//     }
//     return a;
// }

pair<mat,mat> Vs_assigner(int node_x, int node_y, double V_value, mat LHS, mat RHS){

    vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    mat New_LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    mat New_RHS = join_cols(RHS, value);
    return {New_LHS,New_RHS};
}

void print_v(std::vector <mat> const &a) {
   std::cout << "The vector elements are : ";

   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
   std::cout << "The size of the vector is " << a.size();
}
vec arange(double tstart, double h, double vec_size){
    vec time = zeros(vec_size,1);
    time[0] = tstart;
    for(int i=1; i<vec_size; ++i)
    {
        time[i] = tstart + h;
        tstart = tstart + h;
    }
    return time;
}
int main(int argc, const char ** argv){
    double n = 5001;
    double t_start = 0;
    double t_end = 12e-2;
    int t1 = 0;
    double t = t_end - t_start;
    double h = t/(n-1);
    mat LHS(5,5,fill::randu);
    // mat RHS(5,1,fill::randu);
    // LHS.print("The old LHS is:");
    // RHS.print("The old RHS is:");
    // // Assigning a voltage source
    // pair<mat,mat> Matrix = Vs_assigner(3,0,3,LHS,RHS);
    // LHS = Matrix.first;
    // RHS = Matrix.second;
    // // location of voltage value for transient simulation
    // int V_locate1 = RHS.n_cols;
    // // End of assignment
    // LHS.print("The New LHS is:");
    // RHS.print("The New RHS is:");
    
    
    vec test = arange(t_start, h, n);
    test.print("The time vector:");
    return 0;
}
    //Initialize the random generator 
    // arma::arma_rng::set_seed_random();

    // //Create a 4x4 random matrix and print it on the screeen
    // arma::Mat<double> A = arma::randu(4,4);
    // std::cout << "A:\n" << A << "\n";

    // //Multiply A with this transpose:
    // std::cout << "A * A.t() = \n";
    // std::cout << A * A.t() << "\n";
    // //Access/Modify rows and columns from the array:
    // A.row(0) = A.row(1) + A.row(3);
    // A.col(3).zeros();
    // std::cout << "add rows 1 and 3, store result in row 0, also fill 4th column with zeros:\n";
    // std::cout << "A:\n" << A << "\n";

    // //Create a new diagonal matrix using the main diagonal of A:
    // arma::Mat<double>B = arma::diagmat(A);
    // std::cout << "B:\n" << B << "\n";

    // //Save matrices A and B:
    // A.save("A_mat.txt", arma::arma_ascii);
    // B.save("B_mat.txt", arma::arma_ascii);