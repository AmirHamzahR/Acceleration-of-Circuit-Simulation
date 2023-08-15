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
    8) N-MOS Transistor - NMOS_assigner(number, node_vd, node_vg, node_vs, node_vb, h, solution, LHS, RHS, mode);
    9) P-MOS Transistor - PMOS_assigner(number, node_vs, node_vg, node_vd, node_vb, h, solution, LHS, RHS, mode);

    The linear components could be assigned inside the main function while the non-linear and dynamic components can be assigned inside the DynamicNonLinear function.

    V_pulse also has a special assignment where the normal voltage source needs to be assigned with the initial voltage (V1) inside the RHS_locate list. 
    The V_pulse function will then be called inside the transient simulation loop.
*/

/*  TOTAL NUMBER OF NODES EXCLUDING GROUND
    Two port components such as resistors, initially adds 2 nodes. If more than 1 component is added, it then adds 1 node per component.
    Number of MOSFETs and cascaded levels are assigned above too. External nodes are the nodes which  */
int const code = 1; // choosing code for IV_curve
double W = 500e-9;
double L = 50e-9;

int const external_nodes = 4; // Number of external nodes (excluding ground and ring oscillator loop nodes)
int const external_mosfets = 0; // Number of standalone mosfets (excluding mosfets from ring oscillator)
int const cascaded_level = 0; // Number of cascaded ring oscillators 
int const supply_voltage_node = 1; // supply voltage node for the ring oscillator
int const no_of_mosfets = external_mosfets;// Total number of MOSFETs

int const T_nodes = external_nodes + 4*no_of_mosfets;

#include "Transient_code.h"

/*  The line above the code means that the section can be changed by the user which analyses circuit simulation

    --fixed--           : can't be edited for circuit simulation purposes, but can be edited for code debugging.
    --can be changed--  : can be edited for circuit simulation purposes and edited for code debugging.
    
*/

// Assigning the stamp matrices for dynamic and non-linear components
std::pair<arma::mat,arma::mat> DynamicNonLinear(arma::mat &LHS, arma::mat &RHS, arma::mat solution, double h, int mode){
    // (Diode_assigner, PMOS_assigner, NMOS_assigner, C_assigner, RingOscillatorStages)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    NMOS_assigner(1,3,2,0,0,W,L,h,solution,LHS,RHS,mode);
    
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
    
    int i = 0;
    int n = 5001; // Number of iterations 
    // Defining the voltage step settings for VDS
    double V1_start = 0;
    double V1_end = 30;
    double V1_step = 0.1; // Voltage increment
    double h = (V1_end-V1_start)/(n);
    
    // Defining the voltage settings for VGS
    double VGS_value = 10;

    // voltage step vector to be inputted in plot for python analysis
    arma::mat volt_step = arange(V1_start,h,n);
    
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // ASSIGNING THE STAMPS TO THE LHS AND RHS MATRICES

    // default state
    arma::mat LHS = arma::zeros(Maxi,Maxj); // LHS matrix
    arma::mat RHS = arma::zeros(Maxi,1); // RHS matrix

    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE RESISTOR STAMP (R_assigner)
    R_assigner(1,2,220,LHS,RHS);
    R_assigner(3,4,10,LHS,RHS);


    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ASSIGNING THE CURRENT STAMP (Is_assigner, VCCS_assigner)


    /*--------------------------------------------can be changed-------------------------------------------------*/
    // Assigning DC voltage sources
    Vs_assigner(1,0,VGS_value,LHS,RHS); // vgs

    // Assigning the stamps that would affect the RHS in transient simulation 
    // (only for  time-dependent voltage, e.g. pulse voltages)
    std::vector<double> RHS_locate = {
        // Assigning the voltage matrix on LHS and RHS for the pulse voltage
        Vs_assigner(4,0,V1_start,LHS,RHS) // vds
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
    /*----------------------------------------------fixed--------------------------------------------------------*/
    // The solution csv that is going to be plotted which contains the values of nodal voltages
    // and voltage source currents
    arma::vec id_csv = arma::ones(n,1);

    int iter = 0;
    solution = NewtonRaphson_system(init_LHS,init_RHS, LHS, RHS, solution,h,mode);
    id_csv[0] = NMOS_assigner(1,3,2,0,0,W,L,h,solution,LHS,RHS,mode);
    /*--------------------------------------------can be changed-------------------------------------------------*/
    // ADDING VOLTAGE STEPPING LOOP  
    while(i < n){
        
        LHS = init_LHS;
        RHS = init_RHS;
        
        std::vector<double> RHS_value = {
            Vstep_Source(V1_start,V1_step,h) // VDS
        };
        RHS = RHS_update(RHS_locate, init_RHS, RHS_value);
        // Calling the Newton-Raphson system here
        solution = NewtonRaphson_system(init_LHS,init_RHS, LHS, RHS, solution,h,mode);
        
        // Assigning the variables that will be plotted and analysed as seen in a circuit simulator
        id_csv[i] = NMOS_assigner(1,3,2,0,0,W,L,h,solution,LHS,RHS,mode);
        
        i++;
        iter += Maxi;
        
    }
    /*-----------------------------------------------------------------------------------------------------------*/
    // SAVING THE SOLUTION AND TIME MATRICES INTO CSV FILES
    std::ofstream file("id.csv");
    file << "id_matrix" << std::endl;
    id_csv.save(file, arma::csv_ascii);
    file.close();
    std::ofstream file3("VDS_step.csv");
    file3 << "VDS_step" << std::endl;
    volt_step.save(file3, arma::csv_ascii);
    file3.close();

    return 0;
    /*-----------------------------------------------------------------------------------------------------------*/
}