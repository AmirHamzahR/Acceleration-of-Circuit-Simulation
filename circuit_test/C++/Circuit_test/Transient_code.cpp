#include <iostream>
#include "armadillo"
#include <fstream>
#include <armadillo>
#include <tuple>
#include <cmath>
using namespace arma;
using namespace std;

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

// extra -  adding arange() function here to return the time vector
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

// Capacitor stamp assigner
pair<mat,mat> C_assigner(int node_x,int node_y,double C, double h, mat LHS, mat RHS, mat init){
    // this if else statment uses trapezoidal formula
    int M_cols = LHS.n_cols;
    int M_rows = LHS.n_rows;

    double x{0};
    double x1{0};
    
    if(h>0){
        x = C/h;
        if(node_x == 0){
            x1 = C*(init(node_y-1,0))/h;
        }else if(node_y == 0){
            x1 = C*(init(node_x-1,0))/h;
        }else{
            x1 = C*(init(node_x-1,0)-init(node_y-1,0))/h;
        }
    }
    else{
        x = 0;
        x1 = 0;
    }

    // Matrix stamp for a capacitor on RHS
    mat a = Is_assigner(node_x,node_y,x1,M_rows,M_cols);
    // Matrix stamp for a capacitor on LHS
    mat b = R_assigner(node_x,node_y,x,M_rows,M_cols);
    
    return {b, a};
}


// 3) Add pulse voltage system

pair<double,double> V_pulse(double V1, double V2, double t1,double td, double tr, double tf, double tpw, double tper){
    double V_t1 = V1;
    if((t1 >= 0) && (t1 < td)){
        V_t1 = V1;
    }else if((t1 >= td) && (t1 < (td + tr))){
        V_t1 = V1 + (V2-V1)*(t1-td)/tr;
    }else if((t1 >= (td + tr)) && (t1 < (td + tr + tpw))){
        V_t1 = V2;
    }else if((t1 >= (td + tr + tpw)) && (t1 < (td + tr + tpw + tf))){
        V_t1 = V2 + (V1 - V2)*(t1-(td+tr+tpw))/tf;
    }else if((t1 >= (td + tr + tpw + tf)) && (t1 <= (td + tper))){
        V_t1 = V1;
    }else{
        t1 = td;
    }
    return {V_t1, t1};
}
// 2) Add diode system

pair<mat,mat> Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, mat LHS, mat RHS, mat solution){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;
    double x{0};
    double x1{0};
    double val_nodex = solution(node_x-1,0);
    double val_nodey = solution(node_y-1,0);

    if(h>0){
        if(node_x == 0){
            x = (Is/VT)*(exp(-val_nodey/VT)) + cd/h;
            x1 = (x*(-val_nodey)-Is*(exp((-val_nodey)/VT)-1));
        }
        else if(node_y == 0){
            x = (Is/VT)*(exp((val_nodex)/VT)) + cd/h;
            x1 = (x*(val_nodex)-Is*(exp((val_nodex)/VT)-1));
        }
        else{
            x = (Is/VT)*(exp((val_nodex-val_nodey)/VT)) + cd/h;
            x1 = (x*(val_nodex-val_nodey)-Is*(exp((val_nodex-val_nodey)/VT)-1));
        }
    }
    else{
        if(node_x == 0){
            x = (Is/VT)*(exp(-val_nodey/VT));
            x1 = (x*(-val_nodey)-Is*(exp((-val_nodey)/VT)-1));
        }
        else if(node_y == 0){
            x = (Is/VT)*(exp((val_nodex)/VT));
            x1 = (x*(val_nodex)-Is*(exp((val_nodex)/VT)-1));
        }
        else{
            x = (Is/VT)*(exp((val_nodex-val_nodey)/VT));
            x1 = (x*(val_nodex-val_nodey)-Is*(exp((val_nodex-val_nodey)/VT)-1));
        }
    }
    // Matrix stamp for a diode on RHS
    mat a = Is_assigner(node_x,node_y,x1,row_size,col_size);
    // Matrix stamp for a diode on LHS
    mat b = R_assigner(node_x,node_y,x,row_size,col_size);
    // a.print();
    // b.print();
    return {b,a};
}

// Assigning the stamp matrices for dynamic and non-linear components
pair<mat,mat> DynamicNonLinear(mat LHS, mat RHS, mat solution, double h){
    // start of components assignment
    // pair<mat,mat> Diode1 = Diode_assigner(2,3,3e-9,0.05,4e-12,h,LHS,RHS,solution);
    pair<mat,mat> C1 = C_assigner(2,0,0.4e-6,h,LHS,RHS,solution);
    pair<mat,mat> C2 = C_assigner(3,0,0.8e-6,h,LHS,RHS,solution);

    vector<mat> delta_LHS = {
        LHS,
        // Diode1.first,
        C1.first,
        C2.first
    };
    vector<mat> delta_RHS = {
        RHS,
        // Diode1.second,
        C1.second,
        C2.second
    };
    mat J_x = mat_sum(delta_LHS);
    mat F_x = mat_sum(delta_RHS);
    // LHS.print();
    // RHS.print();
    // end of components assignment 
    return {J_x, F_x};
    
}

// 4) Add Newton Raphson system
mat NewtonRaphson_system(mat LHS, mat RHS, mat init, double h){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;
    
    double eps_val = 1e-8;
    mat error = zeros(1,1);
    double error_val = 9e9;
    error.row(0) = error_val;
    int iteration_counter{0};
    mat delta = zeros(row_size,1);
    mat solution = init;

    // Assigning the non-linear and dynamic components
    // pair<mat,mat> Matrices = DynamicNonLinear(LHS,RHS,solution,h);
    pair<mat,mat> Matrices = DynamicNonLinear(LHS,RHS,solution,h);
    mat J_x = Matrices.first;
    mat F_x = Matrices.second;

    while((error(0,0) > eps_val) && (iteration_counter < 50)){
        // The way that Newton Raphson algorithm being solved is a bit different maybe in C++
        delta = solve(J_x,(J_x*solution) - F_x);
        error.row(0) = max(abs(delta));
        solution -= delta;

        pair<mat,mat> Matrices2 = DynamicNonLinear(LHS,RHS,solution,h);
        J_x = Matrices2.first;
        F_x = Matrices2.second;
        
        iteration_counter += 1;
        cout << iteration_counter << " ";
        J_x.print();
        F_x.print();
        // solution.print("The solution matrix:");
    }
    return solution;
}

// Saving to csv function
// void save_csv(mat mat_name){
//     ofstream file(".csv");
//     file << "X,Y" << endl;
//     mat_name.save(file, csv_ascii);
//     file.close();
// }

// Main function for the circuit simulation
int main(int argc, const char ** argv){

    // Total number of nodes excluding ground
    int T_nodes{3};
    // Size of matrix
    int Maxi{T_nodes};
    int Maxj{Maxi};

    // TRANSIENT SIMULATION SETTINGS
    // The amount of iterations for the timestep, the higher the more accurate but uses more computing resources
    int n{5001}; // 5001 seems to be the sweet spot
    int i{0};
    // Defining the Time and Timestep for Transient Simulation
    mat X1 = ones(n,1);
    double t_start = 0;
    double t_end = 12e-2;
    double t = t_end - t_start;
    double h = t/(n-1);
    cout << h << "\n";
    // time vector to be inputted in plot for python analysis
    vec time = arange(t_start,h,n);

    // The nodal voltages that are going to be analysed and plotted
    vec volt1 = ones(n,1);
    vec volt2 = ones(n,1);
    vec volt3 = ones(n,1);
    vec volt4 = ones(n,1);
    vec current = ones(n,1);
    

    // ASSIGNING THE RESISTOR STAMP
    vector<mat> Total_LHS = {
        // Default state - R_assigner(1,0,0,Maxi,Maxj)
        R_assigner(2,1,cond(3e3),Maxi,Maxj),
        R_assigner(3,2,cond(1e3),Maxi,Maxj)
    };

    // ASSIGNING THE CURRENT STAMP
    vector<mat> Total_RHS = {
        // Default state - Is_assigner(1,0,0,Maxi,Maxj)
        Is_assigner(1,0,0,Maxi,Maxj)
    };

    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES
    mat LHS = mat_sum(Total_LHS);
    mat RHS = mat_sum(Total_RHS);
    

    // ASSIGNING A VOLTAGE PULSE
    // Pulse voltage settings
    double V1 = 2;
    double V2 = 6;
    double t1 = 0;
    double td = 1e-3;
    double tr = 0.5e-3;
    double tf = 0.2e-3;
    double tpw = 2e-3;
    double tper = 4e-3;

    pair<mat,mat> voltage1 = Vs_assigner(1,0,V1,LHS,RHS);
    LHS = voltage1.first;
    RHS = voltage1.second;

    // location of voltage value for transient simulation
    int V_locate1 = RHS.n_rows;
    // End of assignment a voltage sources

    // Checking the LHS and RHS matrices
    LHS.print("LHS matrix =");
    RHS.print("RHS matrix =");

    // OPERATING POINT ANALYSIS SYSTEM
    // zero as initial condition
    Maxi = RHS.n_rows;
    Maxj = RHS.n_cols;
    mat init = zeros(Maxi,Maxj);
    // OP analysis used as initial condition for next evaluation
    init = NewtonRaphson_system(LHS,RHS,init,0);

    init.print("The OP analysis of the circuit is: ");

    // ADDING TRANSIENT SIMULATION LOOP
    while(i < n){
        // Loop for the pulse voltage source
        t1 = t1 + h;
        pair<double,double> volt_pulse = V_pulse(V1,V2,t1,td,tr,tf,tpw,tper);
        RHS.row(V_locate1-1).col(0) = volt_pulse.first;
        t1 = volt_pulse.second;
        // cout << i << " " << volt_pulse.first << "\n";

        mat ans = NewtonRaphson_system(LHS,RHS,init,h);

        volt1[i] = ans[0];
        volt2[i] = ans[1];
        volt3[i] = ans[2];
        // volt4[i] = ans[3];
        current[i] = ans[3];
        i++;
    }
    ofstream file("V1.csv");
    file << "V1" << endl;
    volt1.save(file, csv_ascii);
    file.close();
    ofstream file1("V2.csv");
    file1 << "V2" << endl;
    volt2.save(file1, csv_ascii);
    file1.close();
    ofstream file2("V3.csv");
    file2 << "V3" << endl;
    volt3.save(file2, csv_ascii);
    file2.close();
    // ofstream file3("V4.csv");
    // file3 << "V4" << endl;
    // volt4.save(file3, csv_ascii);
    // file3.close();
    ofstream file4("current.csv");
    file4 << "current" << endl;
    current.save(file4, csv_ascii);
    file4.close();
    ofstream file5("time.csv");
    file5 << "time" << endl;
    time.save(file5, csv_ascii);
    file5.close();


    return 0;
}