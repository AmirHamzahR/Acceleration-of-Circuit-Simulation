/*  This code could run a Direct Current Operating Point (DC OP) Analysis and Transient Analysis of a circuit. The circuit could be consisting of
    resistors, capacitors, diodes, current source, voltage source, voltage-controlled current source (VCCS), pulsed voltage source, and N-MOS transistor. 

    The assignments are given as following (with solution, LHS, and RHS being constant):
    1) Resistor - R_assigner(node_x, node_y, cond(R), LHS, RHS);
    2) Capacitors - C_assigner(node_x, node_y, C, h, LHS, RHS, solution, mode);
    3) Diodes - Diode_assigner( node_x, node_y, Is, VT, cd, h, LHS, RHS, solution, mode);
    4) Current Source - Is_assigner( node_x, node_y, I, LHS, RHS);
    5) Voltage Source - Vs_assigner(node_x, node_y, V_value, LHS, RHS);
    6) VCCS - VCCS_assigner(node_x, node_y, node_cx, node_cy, R, LHS);
    7) Pulsed Voltage Source - V_pulse(V1, V2, t1, td, tr, tf, tpw, tper, h);
    8) N-MOS Transistor - fet_assigner(number, node_vd, node_vg, node_vs, node_vb, h, solution, LHS, RHS, mode); (currently maximum of 2 transistors)

    The linear components could be assigned inside the main function while the non-linear and dynamic components can be assigned inside the DynamicNonLinear function.

    V_pulse also has a special assignment where the normal voltage source needs to be assigned with the initial voltage (V1) inside the RHS_locate list. 
    The V_pulse function will then be called inside the transient simulation loop.
*/

/*  TOTAL NUMBER OF NODES EXCLUDING GROUND
    Two port components such as resistors, initially adds 2 nodes. If more than 1 component is added, it then adds 1 node per component.
    Each NMOS added, adds 4 internal nodes. The other nodes can be added separately in the same equation. */
int const T_nodes = 3 /*External nodes*/ + 4*0 /*NMOS (internal nodes)*/ + 4*1 /*PMOS (internal nodes)*/;

#include "Transient_code.h"

/*  The line above the code means that the section can be changed by the user which analyses circuit simulation

    --fixed--           : can't be edited for circuit simulation purposes, but can be edited for code debugging.
    --can be changed--  : can be edited for circuit simulation purposes and edited for code debugging.
    
*/

// Assigning the stamp matrices for dynamic and non-linear components
std::pair<arma::mat,arma::mat> DynamicNonLinear(arma::mat LHS, arma::mat RHS, arma::mat solution, double h, int mode){
    // (Diode_assigner, fet_assigner, C_assigner)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    PMOS_assigner(1, 2, 1, 0, 0, h, solution, LHS, RHS, mode);
    arma::mat J_x = LHS;
    arma::mat F_x = RHS;
    
    return {J_x,F_x};
}

// Main function for the circuit simulation
int main(int argc, const char ** argv){
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Size of matrix
    int Maxi{T_nodes}; // defined in header file (Transient_code.h)
    int Maxj{Maxi};

    /*--------------------------------------------can be changed-------------------------------------------------*/

    // TRANSIENT SIMULATION SETTINGS
    // The amount of iterations for the timestep, the higher the more accurate but uses more computing resources
    int n = 10001; // 5001 seems to be the sweet spot
    int i = 0;
    // Defining the Time and Timestep for Transient Simulation
    arma::mat X1 = arma::ones(n,1);
    double t_start = 0;
    double t_end = 20e-3;
    double t = t_end - t_start;
    double h = t/(n-1);
    std::cout << h << "\n";
    // time vector to be inputted in plot for python analysis
    arma::mat time = arange(t_start,h,n);
    
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    // default state
    arma::mat LHS = arma::zeros(Maxi,Maxj); // LHS matrix
    arma::mat RHS = arma::zeros(Maxi,1); // RHS matrix

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE RESISTOR STAMP (R_assigner)
    R_assigner(3,2,1e3,LHS,RHS);

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE CURRENT STAMP (Is_assigner, VCCS_assigner)



    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE VOLTAGE SOURCES (Vs_assigner)
    // Pulse voltage settings
    double t1 = 0; // time used for the loop
    double V1 = 0;
    double V2 = -5;
    double td = 1e-3;
    double tr = 1e-3;
    double tf = 1e-3;
    double tpw = 0;
    double tper = 2e-3;

    // Assigning DC voltage sources
    // Vs_assigner(1,0,2,LHS,RHS); // VDD
    Vs_assigner(3,0,-0.5,LHS,RHS);

    // Assigning the stamps that would affect the RHS in transient simulation 
    // (only for  time-dependent voltage, e.g. pulse voltages)
    std::vector<double> RHS_locate = {
        // Assigning the voltage matrix on LHS and RHS for the pulse voltage
        Vs_assigner(1,0,V1,LHS,RHS)
    };
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // Checking the LHS and RHS matrices
    LHS.print("LHS matrix =");
    RHS.print("RHS matrix =");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // OPERATING POINT ANALYSIS SYSTEM
    int mode = 0; // 0 to do OP analysis, 1 to do transient simulation
    // zero as initial condition
    Maxi = RHS.n_rows;
    Maxj = RHS.n_cols;

    // The initial LHS and RHS values to be used in the NR-algorithm
    arma::mat const init_LHS = LHS;
    arma::mat const init_RHS = RHS;
    
    // OP analysis used as initial condition for next evaluation
    arma::vec solution= arma::zeros(Maxi,Maxj);
    solution = NewtonRaphson_system(init_LHS,init_RHS,LHS,RHS,solution,h,mode);
    solution.print("The OP analysis of the circuit is: ");
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // The solution csv that is going to be plotted which contains the values of nodal voltages
    // and voltage source currents
    arma::vec solution_csv = arma::ones(n*Maxi,1);
    arma::vec Max_I = arma::zeros(1,1);
    Max_I.row(0).col(0) = Maxi;
    int iter = 0;
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ADDING TRANSIENT SIMULATION LOOP (includes V_pulse or any time dependent sources)
    mode = 1;
    while(i < n){
        
        LHS = init_LHS;
        RHS = init_RHS;
        std::vector<double> RHS_value = {
            V_pulse(V1,V2,t1,td,tr,tf,tpw,tper,h)
        };
        RHS = RHS_update(RHS_locate, init_RHS, RHS_value);
        // Calling the Newton-Raphson system here
        solution = NewtonRaphson_system(init_LHS,init_RHS, LHS, RHS, solution,h,mode);
        // Assigning the variables that will be plotted and analysed as seen in a circuit simulator
        for(int a = 0; a<Maxi; a++){
            solution_csv[iter+a] = solution[a];
        }
        i++;
        iter += Maxi;
    }
    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES
    std::ofstream file("solution.csv");
    file << "X_matrix" << std::endl;
    solution_csv.save(file, arma::csv_ascii);
    file.close();
    std::ofstream file2("MaxI.csv");
    file2 << "Max_I" << std::endl;
    Max_I.save(file2, arma::csv_ascii);
    file2.close();
    std::ofstream file3("time.csv");
    file3 << "time" << std::endl;
    time.save(file3, arma::csv_ascii);
    file3.close();

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}