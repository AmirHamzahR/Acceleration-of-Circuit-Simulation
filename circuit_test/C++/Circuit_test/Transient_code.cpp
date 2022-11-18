#include <iostream>
#include "armadillo"
#include <fstream>
#include <armadillo>
#include <tuple>
#include <cmath>

// Matrix stamps assigner using Modified Nodal Analysis
void DynamicNonLinear(arma::mat &LHS, arma::mat &RHS, arma::mat solution, double h, int mode);
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS);
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS);
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
void C_assigner(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);


// Sum the matrices inside the vector
arma::mat mat_sum(std::vector<arma::mat> vector_of_matrices){
    arma::mat a = vector_of_matrices[0];
    for (long unsigned int i = 1; i < vector_of_matrices.size(); ++i) {
        a = a + vector_of_matrices[i];
    }
    return a;
}

// Testing LU solve using the triangular pivoting system
arma::mat LU_solve(arma::mat A, arma::mat b){
    arma::mat L, U, P;
    arma::lu(L,U,P,A);
    arma::vec x1 = arma::solve(trimatu(U), solve(trimatl(L), P*b));

    return x1;
}

// Changing resistance to conductance
double cond(double R){
    return 1/R;
}

// extra -  adding arange() function here to return the time vector
arma::vec arange(double tstart, double h, double vec_size){
    arma::vec time = arma::zeros(vec_size,1);
    for(int i=0; i<vec_size; ++i)
    {
        time[i] = tstart;
        tstart = tstart + h;
    }
    return time;
}

// create resistor matrix stamp
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    arma::mat a = arma::zeros(maxi,maxj);
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
    LHS = LHS + a;
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = 1;
    arma::mat a = arma::zeros(maxi,maxj);
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
    RHS = RHS + a;
}
// branch extender function
arma::mat branch_ext(arma::mat M, int node_x, int node_y){
    int M_cols = M.n_cols;
    int M_rows = M.n_rows;
    arma::mat va = arma::zeros(M_rows,1);
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

    arma::mat zero_ext = arma::zeros(1,1);
    arma::mat ha = va.as_row();
    arma::mat haz = arma::join_rows(ha,zero_ext);

    arma::mat M1 = arma::join_rows(M,va);
    arma::mat M2 = arma::join_cols(M1,haz);

    return M2;
}

// Voltage source stamp assigner
void Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);
}

// Capacitor stamp assigner
void C_assigner(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode){
    // this if else statment uses trapezoidal formula
    double x = 0;
    double x1 = 0;

    if(mode>0){
        x = C/h;
        if(node_x == 0){
            x1 = C*(solution(node_y-1,0))/h;
        }else if(node_y == 0){
            x1 = C*(solution(node_x-1,0))/h;
        }else{
            x1 = C*(solution(node_x-1,0)-solution(node_y-1,0))/h;
        }
    }
    else{
        x = 0;
        x1 = 0;
    }
    // Matrix stamp for a capacitor on LHS
    R_assigner(node_x,node_y,x,LHS,RHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x,node_y,x1,LHS,RHS);
}

// Voltage pulse assigner
std::pair<double,double> V_pulse(double V1, double V2, double t1,double td, double tr, double tf, double tpw, double tper){
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

// Diode stamp assigner
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;

    double x = 0;
    double x1 = 0;
    double x2 = 0;
    double val_nodex = solution(node_x-1,0);
    double val_nodey = solution(node_y-1,0);

    if(mode>0){
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
    Is_assigner(node_x,node_y,x1,LHS,RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x,node_y,x,LHS,RHS);
}

// Newton Raphson system solver for non-linear and dynamic elements
arma::mat NewtonRaphson_system(arma::mat LHS, arma::mat RHS, arma::mat &init, double h, int mode){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;
    
    double eps_val = 1e-8;
    arma::mat error = arma::zeros(1,1);
    double error_val = 9e9;
    error.row(0) = error_val;
    int iteration_counter= 0;
    arma::mat delta = arma::zeros(row_size,1);
    arma::mat solution = init;
    
    while((error(0,0) > eps_val) && (iteration_counter < 5)){

        DynamicNonLinear(LHS,RHS,solution,h,mode);
        delta = LU_solve(LHS,(LHS*solution) - RHS);
        error.row(0) = arma::max(arma::abs(delta));
        solution -= delta;
        
        iteration_counter += 1;
        // std::cout << iteration_counter;
    }
    return solution;
}

// Assigning the stamp matrices for dynamic and non-linear components
void DynamicNonLinear(arma::mat &LHS, arma::mat &RHS, arma::mat solution, double h, int mode){
    Diode_assigner(2,3,2.7e-9,0.05,4e-12,h,LHS,RHS,solution,mode);
    C_assigner(3,0,0.4e-6,h,LHS,RHS,solution,mode);
    C_assigner(4,0,0.8e-6,h,LHS,RHS,solution,mode);
}

// Main function for the circuit simulation
int main(int argc, const char ** argv){

    // Total number of nodes excluding ground
    int T_nodes = 4;
    // Size of matrix
    int Maxi{T_nodes};
    int Maxj{Maxi};

    // TRANSIENT SIMULATION SETTINGS
    // The amount of iterations for the timestep, the higher the more accurate but uses more computing resources
    int n = 5001; // 5001 seems to be the sweet spot
    int i = 0;
    // Defining the Time and Timestep for Transient Simulation
    arma::mat X1 = arma::ones(n,1);
    double t_start = 0;
    double t_end = 12e-2;
    double t = t_end - t_start;
    double h = t/(n-1);
    std::cout << h << "\n";
    // time vector to be inputted in plot for python analysis
    arma::mat time = arange(t_start,h,n);

    // The nodal voltages that are going to be analysed and plotted
    arma::vec volt1 = arma::ones(n,1);
    arma::vec volt2 = arma::ones(n,1);
    arma::vec volt3 = arma::ones(n,1);
    arma::vec volt4 = arma::ones(n,1);
    arma::vec current = arma::ones(n,1);
    
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES
    // default state
    arma::mat LHS = arma::zeros(Maxi,Maxj);
    arma::mat RHS = arma::zeros(Maxi,1);

    // ASSIGNING THE RESISTOR STAMP
    R_assigner(2,1,cond(3e3),LHS,RHS);
    R_assigner(4,3,cond(1e3),LHS,RHS);

    // ASSIGNING THE CURRENT STAMP

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
    
    // Assigning the voltage matrix on LHS and RHS
    Vs_assigner(1,0,V1,LHS,RHS);

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
    arma::mat init = arma::zeros(Maxi,Maxj);
    // OP analysis used as initial condition for next evaluation
    int mode = 0; // 0 to do OP analysis
    init = NewtonRaphson_system(LHS,RHS,init,h,mode);

    init.print("The OP analysis of the circuit is: ");
    // ADDING TRANSIENT SIMULATION LOOP
    while(i < n){
        // Loop for the pulse voltage source
        t1 = t1 + h;
        std::pair<double,double> volt_pulse = V_pulse(V1,V2,t1,td,tr,tf,tpw,tper);
        RHS.row(V_locate1-1).col(0) = volt_pulse.first;
        t1 = volt_pulse.second;
        // Calling the Newton-Raphson system here
        mode = 1; // 1 to do transient simulation
        init = NewtonRaphson_system(LHS,RHS,init,h,mode);
        // Assigning the variables that will be plotted and analysed as seen in a circuit simulator
        volt1[i] = init[0];
        volt2[i] = init[1];
        volt3[i] = init[2];
        volt4[i] = init[3];
        current[i] = init[4];
        i++;
    }
    std::ofstream file("V1.csv");
    file << "V1" << std::endl;
    volt1.save(file, arma::csv_ascii);
    file.close();
    std::ofstream file1("V2.csv");
    file1 << "V2" << std::endl;
    volt2.save(file1, arma::csv_ascii);
    file1.close();
    std::ofstream file2("V3.csv");
    file2 << "V3" << std::endl;
    volt3.save(file2, arma::csv_ascii);
    file2.close();
    std::ofstream file3("V4.csv");
    file3 << "V4" << std::endl;
    volt4.save(file3, arma::csv_ascii);
    file3.close();
    std::ofstream file4("current.csv");
    file4 << "current" << std::endl;
    current.save(file4, arma::csv_ascii);
    file4.close();
    std::ofstream file5("time.csv");
    file5 << "time" << std::endl;
    time.save(file5, arma::csv_ascii);
    file5.close();


    return 0;
}