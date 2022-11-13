#include <iostream>
#include "armadillo"
#include <armadillo>
#include <tuple>
using namespace arma;
using namespace std;

// create resistor matrix stamp
mat R_assigner(double node_x, double node_y, double R, double maxi, double maxj){
    mat a = zeros(maxi,maxj);
    if(node_x == 0){
        a.row(node_y-1).col(node_y-1) = R;
    }
    else if(node_y == 0){
        a.row(node_x-1).col(node_x-1) = R;
    }
    else{
        a.row(node_x-1).col(node_x-1) = R;
        a.row(node_x-1).col(node_y-1) = -R;
        a.row(node_y-1).col(node_x-1) = -R;
        a.row(node_y-1).col(node_y-1) = R;
    }
    return a;
}

// create a current matrix stamp
mat Is_assigner(double node_x, double node_y, double I, int maxi, int maxj){
    maxj = 1;
    mat a = zeros(maxi,maxj);
    if(node_x == 0){
        a.row(node_y-1).col(0) = -I;
    }
    else if(node_y == 0){
        a.row(node_x-1).col(0) = I;
    }
    else{
        a.row(node_x-1).col(0) = I;
        a.row(node_y-1).col(0) = -I;
    }
    return a;
}
// branch extender function
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

// Voltage source stamp assigner
pair<mat,mat> Vs_assigner(int node_x, int node_y, double V_value, mat LHS, mat RHS){

    vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    mat New_LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    mat New_RHS = join_cols(RHS, value);
    return {New_LHS,New_RHS};
}

mat mat_sum(vector<mat> vector_of_matrices){
    mat a = vector_of_matrices[0];
    for (int i = 1; i < vector_of_matrices.size(); ++i) {
        a = a + vector_of_matrices[i];
    }
    return a;
}

mat LU_solve(mat A, mat b){
    mat L, U, P;
    lu(L,U,P,A);
    vec x1 = solve(trimatu(U), solve(trimatl(L), P*b));

    return x1;
}

double cond(double R){
    return 1/R;
}

int main(int argc, const char ** argv){
    // Total number of nodes
    int T_nodes{5};
    // Size of matrix
    int Maxi{T_nodes};
    int Maxj{Maxi};

    // assigning the resistor stamps
    vector<mat> Total_LHS = {
        // Default state - R_assigner(0,0,0,Maxi,Maxj)
        R_assigner(2,1,cond(2),Maxi,Maxj),
        R_assigner(2,0,cond(3),Maxi,Maxj),
        R_assigner(3,2,cond(4),Maxi,Maxj),
        R_assigner(3,0,cond(2),Maxi,Maxj),
        R_assigner(3,4,cond(1),Maxi,Maxj)
    };

    // assigning the current stamps
    vector<mat> Total_RHS = {
        // Default state - Is_assigner(0,0,0,Maxi,Maxj)
        Is_assigner(4,0,3,Maxi,Maxj)
    };

    mat LHS = mat_sum(Total_LHS);
    LHS.print("Resistor matrix =");
    mat RHS = mat_sum(Total_RHS);
    RHS.print("Current matrix =");

    // Assigning a voltage source
    pair<mat,mat> Matrix = Vs_assigner(1,0,3,LHS,RHS);
    LHS = Matrix.first;
    RHS = Matrix.second;
    // // location of voltage value for transient simulation
    int V_locate1 = RHS.n_cols;
    // End of assignment a voltage sources

    mat x = solve(LHS,RHS);

    x.print("The OP analysis of the circuit is: ");

    return 0;
}