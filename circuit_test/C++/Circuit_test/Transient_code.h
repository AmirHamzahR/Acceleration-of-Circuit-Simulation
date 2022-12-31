#ifndef Transient_code
#define Transient_code
#include <iostream>
#include "armadillo"
#include <fstream>
#include <armadillo>
#include <tuple>
#include <cmath>

// Matrix stamps assigner using Modified Nodal Analysis
std::pair<arma::mat,arma::mat> DynamicNonLinear(arma::mat LHS, arma::mat RHS, arma::mat solution, double h, int mode);
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS);
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS);
double Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS);
void C_assigner(int node_x,int node_y,double C, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode);
void VCCS_assigner(int node_x,int node_y,int node_cx,int node_cy,double R,arma::mat &LHS);
void NMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double h, arma::mat &solution, arma::mat &LHS, arma::mat &RHS,  int mode);

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

// assigning the matrix stamps for the VCCS
void VCCS_assigner(int node_x,int node_y,int node_cx,int node_cy,double R,arma::mat &LHS){
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    arma::mat a = arma::zeros(maxi,maxj);
    
    if(node_x == 0){
        if(node_cx == 0){
            if(node_cy>0){
            a.row(node_y-1).col(node_cy-1) = R;
            }
        }
        else if(node_cy == 0){
            if(node_cx > 0){ 
            a.row(node_y-1).col(node_cx-1) = -R;
            }
        }
        else{ 
        a.row(node_y-1).col(node_cx-1) = -R;
        a.row(node_y-1).col(node_cy-1) = R;
        }
    }
    else if(node_y == 0){
        if(node_cx == 0){
            if(node_cy > 0){
                a.row(node_x-1).col(node_cy-1) = -R;
            }
        }
        else if(node_cy == 0){
            if(node_cx > 0){
                a.row(node_x-1).col(node_cx-1) = R;
            }
        }
        else{
            a.row(node_x-1).col(node_cx-1) = R;
            a.row(node_x-1).col(node_cy-1) = -R;
        }
    }
    else{
        if(node_cx == 0){
            if(node_cy > 0){
                a.row(node_x-1).col(node_cy-1) = -R;
                a.row(node_y-1).col(node_cy-1) = R;
            }
        }
        else if(node_cy == 0){
            if(node_cx > 0){
                a.row(node_x-1).col(node_cx-1) = R;
                a.row(node_y-1).col(node_cx-1) = -R;
            }
        }
        else{
            a.row(node_x-1).col(node_cx-1) = R;
            a.row(node_x-1).col(node_cy-1) = -R;
            a.row(node_y-1).col(node_cx-1) = -R;
            a.row(node_y-1).col(node_cy-1) = R;
        }
    }
    LHS = LHS + a;
}

// create resistor matrix stamp
void R_assigner(double node_x, double node_y, double R, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = LHS.n_rows;
    double x = 0;
    arma::mat a = arma::zeros(maxi,maxj);

    if(R == 0)
        x = 0;
    else
        x = cond(R);

    if((node_x == 0) && (node_y == 0)){
        x = 0;
    }
    else{
        if(node_x == 0){
            a.row(node_y-1).col(node_y-1) = x;
        }
        else if(node_y == 0){
            a.row(node_x-1).col(node_x-1) = x;
        }
        else{
            a.row(node_x-1).col(node_x-1) = x;
            a.row(node_x-1).col(node_y-1) = -x;
            a.row(node_y-1).col(node_x-1) = -x;
            a.row(node_y-1).col(node_y-1) = x;
        }
    }
    LHS = LHS + a;
}

// create a current matrix stamp
void Is_assigner(double node_x, double node_y, double I, arma::mat &LHS, arma::mat &RHS){
    int maxi = LHS.n_cols;
    int maxj = 1;
    arma::mat a = arma::zeros(maxi,maxj);
    if((node_x == 0) && (node_y == 0)){
        I = 0;
    }
    else{
        if(node_x == 0){
            a.row(node_y-1).col(0) = I;
        }
        else if(node_y == 0){
            a.row(node_x-1).col(0) = -I;
        }
        else{
            a.row(node_x-1).col(0) = -I;
            a.row(node_y-1).col(0) = I;
        }
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

// Inductor model - not yet complete, quite hard to add in the circuit 
// double Ind_assigner(int node_x, int node_y, double L, double h, double RHS_locate, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode){
//     double Ind = 0;
//     if(mode > 0){
//         Ind = -L/h*solution(RHS_locate-1,0);
//         LHS.row(LHS.n_rows-1).col(LHS.n_cols-1) = -L/h;
//     }
//     else{
//     arma::mat a = arma::zeros(1,1);
//         // Extending the branch at the LHS matrix
//         LHS = branch_ext(LHS, node_x, node_y);
//         LHS.row(LHS.n_rows-1).col(LHS.n_cols-1) = 0;
//         // Assigning the value at RHS
//         RHS = arma::join_cols(RHS, a);
//         // location of inductor value for transient simulation
        
//         Ind = RHS.n_rows;
//     }
//     return Ind;
// }

void NMOS_assigner(int number, int node_vd, int node_vg, int node_vs, int node_vb, double h, arma::mat &solution, arma::mat &LHS, arma::mat &RHS,  int mode){
    // # the position of the drain, gate, source, and base voltages are hard-coded for the transistor model
    // # we are using a discrete MOSFET, so the vb is connected to the source terminal

    double vd = 0;
    double vg = 0;
    double vs = 0;
    double vb = 0;

    // # The settings for the large signal analysis model
    // # uses node_vd as starting reference node for the simulation
    if(node_vd == 0){
        vd = 0;
    }else
        vd = solution(node_vd-1,0);
        
    if(node_vg == 0){
        vg = 0;
    }else
        vg = solution(node_vg-1,0);
    
    if(node_vs == 0){
        vs = 0;
    }else
        vs = solution(node_vs-1,0);
        
    if(node_vb == 0){
        vb = 0;
    }else
        vb = solution(node_vb-1,0);
    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // depletion mode voltages
    double vbd = vb - vd;
    double vgd = vg - vd;


    double id = 0;
    double gds = 0;
    double gm = 0;
    double gmb = 0;
    double W = 100e-6; // default is 100um
    double L = 100e-6; // default is 100um
    double Ld = 0;
    double Leff = L-2*Ld;
    double kp = 2e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0;
    
    double Beta = (mCox)*(W/Leff);
    double gamma = 0;
    double phi = 0.65;
    double vt0 = 0; // default is 0
    double vt = 0;
    double I_DSeq = 0;

    double CGSO = 4e-11;
    double CGDO = 4e-11;
    double CGBO = 2e-10;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * W;
    double CBD = 6e-17; // typical value for CBD
    double CBS = 6e-17; // typical value for CBS

    // # the settings for fet model based on the large signal analysis
    R_assigner(node_vd,T_nodes-(4*number)+2,1,LHS,RHS); // # RD
    R_assigner(node_vg,T_nodes-(4*number)+1,1,LHS,RHS); // # RG
    R_assigner(T_nodes-(4*number)+4,node_vs,1,LHS,RHS); // # RS
    R_assigner(T_nodes-(4*number)+3,node_vb,1,LHS,RHS); // # RB
    Diode_assigner(T_nodes-(4*number)+3,T_nodes-(4*number)+2,1e-14,0.05,CBD,h,LHS,RHS,solution,mode); // # Diode BD
    Diode_assigner(T_nodes-(4*number)+3,T_nodes-(4*number)+4,1e-14,0.05,CBS,h,LHS,RHS,solution,mode); // # Diode BS
    
    // For the conductances
    // if(vds>=0){
        vt = vt0 + gamma*((sqrt(phi-vbs)-sqrt(phi))); // # already taking into account the body effect of MOSFETs
    // }else{
    //     vt = vt0 + gamma*((sqrt(phi-vbd)-sqrt(phi))); // # already taking into account the body effect of MOSFETs
    // }
    // gm = Beta*vds;
    // gds = Beta*(vgs-vt-vds);
    // if(vds>=0){
        if ((vds <= (vgs-vt)) && (vgs >= vt)){ // # the transistor is in linear
            id = Beta*(vgs-vt-vds/2)*vds*(1+LAMBDA*vds);
            gds = Beta*(1+LAMBDA*vds)*(vgs-vt-vds)+Beta*LAMBDA*vds*(vgs-vt-vds/2);
            gm = Beta*(1+LAMBDA*vds)*vds;
            gmb = gm*gamma/(2*sqrt(phi-vbs));
        }else if ((vds >= (vgs-vt)) && (vgs >= vt)){ // # the transistor is in saturation
            id = (Beta/2)*pow((vgs - vt),2) * (1+LAMBDA*vds);
            gds = (Beta/2)*LAMBDA*pow((vgs-vt),2);
            gm = Beta*(1+LAMBDA*vds)*(vgs-vt);
            gmb = gm*gamma/(2*sqrt(phi-vbs));
        }else{ // # the transistor is in cutoff
            id = 0;
            gds = 0;
            gm = 0;
            gmb = 0;
        }
    // }else{ // For depletion mode
    //     vbs = vbd;
    //     vgs = vgd;
    //     vds = -vds;
    //     if ((vds < (vgs-vt))){ // # the transistor is in linear
    //         id = Beta*(vgs-vt-1/2*vds)*vds*(1+LAMBDA*vds);
    //         gds = Beta*(1+LAMBDA*vds)*(vgs-vt-vds)+Beta*LAMBDA*vds*(vgs-vt-1/2*vds);
    //         gm = Beta*(1+LAMBDA*vds)*vds;
    //         gmb = gm*gamma/(2*sqrt(phi-vbs));
    //     }else if ((vds > (vgs-vt)) && ((vgs-vt)>0)){ // # the transistor is in saturation
    //         id = (Beta/2)*pow((vgs - vt),2) * (1+LAMBDA*vds);
    //         gds = (Beta/2)*LAMBDA*pow((vgs-vt),2);
    //         gm = Beta*(1+LAMBDA*vds)*(vgs-vt);
    //         gmb = gm*gamma/(2*sqrt(phi-vbs));
    //     }else if((vgs-vt)<0){ // # the transistor is in cutoff
    //         id = 0;
    //         gds = 0;
    //         gm = 0;
    //         gmb = 0;
    //     }
    //     id = -id;
    // }

    // might need to change the capacitor values according to the textbook model which takes in the linear, saturation, and cutoff regions
    C_assigner(T_nodes-(4*number)+1,T_nodes-(4*number)+4,CGS,h,LHS,RHS,solution,mode); // # Capacitor GSO
    C_assigner(T_nodes-(4*number)+1,T_nodes-(4*number)+3,CGB,h,LHS,RHS,solution,mode); // # Capacitor GBO
    C_assigner(T_nodes-(4*number)+4,T_nodes-(4*number)+3,CBS,h,LHS,RHS,solution,mode); // # Capacitor BSO

    I_DSeq = id - gds*vds - gm*vgs - gmb*vbs/*/*  /**/; // # 10.190 equation
    // Is_assigner(node_vs,node_vd,id,LHS,RHS);
    Is_assigner(node_vd,node_vs,I_DSeq,LHS,RHS);
    VCCS_assigner(node_vd,node_vs,node_vb,node_vs,gmb,LHS); // assigning gmb
    R_assigner(node_vd,node_vs,cond(gds),LHS,RHS); // # assigning gds
    VCCS_assigner(node_vd,node_vs,node_vg,node_vs,gm,LHS); // # assigning gm
}

void PMOS_assigner(int number, int node_vs, int node_vg, int node_vd, int node_vb, double h, arma::mat &solution, arma::mat &LHS, arma::mat &RHS,  int mode){
    // # the position of the drain, gate, source, and base voltages are hard-coded for the transistor model
    // # we are using a discrete MOSFET, so the vb is connected to the source terminal

    double vd = 0;
    double vg = 0;
    double vs = 0;
    double vb = 0;

    // # The settings for the large signal analysis model
    // # uses node_vd as starting reference node for the simulation
    if(node_vd == 0){
        vd = 0;
    }else
        vd = solution(node_vd-1,0);
        
    if(node_vg == 0){
        vg = 0;
    }else
        vg = solution(node_vg-1,0);
    
    if(node_vs == 0){
        vs = 0;
    }else
        vs = solution(node_vs-1,0);
        
    if(node_vb == 0){
        vb = 0;
    }else
        vb = solution(node_vb-1,0);

    // enhancement mode voltages
    double vgs = vg - vs;
    double vds = vd - vs;
    double vbs = vb - vs;
    // enhancement mode voltages
    double vsg = vs - vg;
    double vsd = vs - vd;
    double vsb = vs - vb;
    // depletion mode voltages
    double vdb = vd - vb;
    double vdg = vd - vg;


    double id = 0;
    double gds = 0;
    double gm = 0;
    double gmb = 0;
    double W = 100e-6; // default is 100um
    double L = 100e-6; // default is 100um
    double Ld = 0;
    double Leff = L-2*Ld;
    double kp = 2e-5; // default is 2e-5
    double mCox = kp;
    double LAMBDA = 0;
    
    double Beta = (mCox)*(W/Leff);
    double gamma = 0;
    double phi = 0.65;
    double vt0 = 0; // default is 0
    double vt = 0;
    double I_DSeq = 0;

    double CGSO = 4e-11;
    double CGDO = 4e-11;
    double CGBO = 2e-10;
    double CGS = CGSO * W;
    double CGD = CGDO * W;
    double CGB = CGBO * W;
    double CBD = 6e-17; // 6e-17 typical value for CBD
    double CBS = 6e-17; // 6e-17 typical value for CBS

    // # the settings for fet model based on the large signal analysis
    R_assigner(node_vd,T_nodes-(4*number)+2,1,LHS,RHS); // # RD
    R_assigner(node_vg,T_nodes-(4*number)+1,1,LHS,RHS); // # RG
    R_assigner(T_nodes-(4*number)+4,node_vs,1,LHS,RHS); // # RS
    R_assigner(T_nodes-(4*number)+3,node_vb,1,LHS,RHS); // # RB
    Diode_assigner(T_nodes-(4*number)+2,T_nodes-(4*number)+3,1e-14,0.05,CBD,h,LHS,RHS,solution,mode); // # Diode BD
    Diode_assigner(T_nodes-(4*number)+4,T_nodes-(4*number)+3,1e-14,0.05,CGD,h,LHS,RHS,solution,mode); // # Diode BS
    
    // If vds>=0, then the transistor is in enhancement mode, otherwise it is in depletion mode
    // if(vsd>=0){
        vt = vt0 - gamma*((sqrt(phi-vsb)-sqrt(phi))); // # already taking into account the body effect of MOSFETs
    // }else{
    //     vt = vt0 + gamma*((sqrt(phi-vdb)-sqrt(phi))); // # already taking into account the body effect of MOSFETs
    // }
    
    double n_vt = std::abs(vt);
    // if (vsd >= 0){
            if ((vds >= (vgs-vt)) && (vgs <= vt)){ // # the transistor is in linear
                id = Beta*(vsg-n_vt-vsd/2)*vsd*(1+LAMBDA*vsd);
                gds = Beta*(1+LAMBDA*vsd)*(vsg-n_vt-vsd)+Beta*LAMBDA*vsd*(vsg-n_vt-vsd/2);
                gm = Beta*(1+LAMBDA*vsd)*vsd;
                gmb = gm*gamma/(2*sqrt(phi-vsb));
                // std::cout<<"Linear"<<std::endl;
             }else if ((vds <= (vgs-vt)) && (vgs <= vt)){ // # the transistor is in saturation
                id = (Beta/2)*pow((vsg - n_vt),2) * (1+LAMBDA*vsd);
                gds = (Beta/2)*LAMBDA*pow((vsg-n_vt),2);
                gm = Beta*(1+LAMBDA*vsd)*(vsg-n_vt);
                gmb = gm*gamma/(2*sqrt(phi-vsb));
                // std::cout<<"Saturation"<<std::endl;
            }else{ // # the transistor is in cutoff
                id = 0;
                gds = 0;
                gm = 0;
                gmb = 0;
                // std::cout<<"Cutoff"<<std::endl;
            } 
    // }else{ // For depletion mode
    //     vsb = vdb;
    //     vsg = vdg;
    //     vsd = -vsd;
    //     if ((vsd <= (vsg-n_vt)) && (vsg > n_vt)){ // # the transistor is in linear
    //             id = Beta*(vsg-n_vt-vsd/2)*vsd*(1+LAMBDA*vsd);
    //             gds = Beta*(1+LAMBDA*vsd)*(vsg-n_vt-vsd)+Beta*LAMBDA*vsd*(vsg-n_vt-vsd/2);
    //             gm = Beta*(1+LAMBDA*vsd)*vsd;
    //             gmb = gm*gamma/(2*sqrt(phi-vsb));
    //         }else if ((vsd > (vsg-n_vt)) && (vsg > n_vt)){ // # the transistor is in saturation
    //             id = (Beta/2)*pow((vsg - n_vt),2) * (1+LAMBDA*vsd);
    //             gds = (Beta/2)*LAMBDA*pow((vsg-n_vt),2);
    //             gm = Beta*(1+LAMBDA*vsd)*(vsg-n_vt);
    //             gmb = gm*gamma/(2*sqrt(phi-vsb));
    //         }else{ // # the transistor is in cutoff
    //             id = 0;
    //             gds = 0;
    //             gm = 0;
    //             gmb = 0;
    //     }
    //     id = -id;
    // }

    // might need to change the capacitor values according to the textbook model which takes in the linear, saturation, and cutoff regions
    C_assigner(T_nodes-(4*number)+1,T_nodes-(4*number)+4,CGS,h,LHS,RHS,solution,mode); // # Capacitor GSO
    C_assigner(T_nodes-(4*number)+1,T_nodes-(4*number)+3,CGB,h,LHS,RHS,solution,mode); // # Capacitor GBO
    C_assigner(T_nodes-(4*number)+4,T_nodes-(4*number)+3,CBS,h,LHS,RHS,solution,mode); // # Capacitor BSO

    I_DSeq = id - gds*vsd - gm*vsg - gmb*vsb;

    Is_assigner(node_vs,node_vd,I_DSeq,LHS,RHS);
    VCCS_assigner(node_vs,node_vd,node_vs,node_vb,gmb,LHS); // assigning gmb
    R_assigner(node_vd,node_vs,cond(gds),LHS,RHS); // # assigning gds
    VCCS_assigner(node_vs,node_vd,node_vs,node_vg,gm,LHS); // # assigning gm
}

void RingOscillatorStages(arma::mat &LHS, arma::mat &RHS, arma::mat solution, double h, int mode){
    // (Diode_assigner, PMOS_assigner, NMOS_assigner, C_assigner)
    /*--------------------------------------------can be changed-------------------------------------------------*/
    int even = 0;
    int odd = 1;
    for(int i = 1; i <= cascaded_level; i++){
        even += 2;
        PMOS_assigner(odd, supply_voltage_node, i+1, i+2, supply_voltage_node, h, solution, LHS, RHS, mode);
        NMOS_assigner(even, i+2, i+1, 0, 0, h, solution, LHS, RHS, mode);
        odd += 2;
    }
    // PMOS_assigner(1, 1, 2, 3, 1, h, solution, LHS, RHS, mode);
    // NMOS_assigner(2, 3, 2, 0, 0, h, solution, LHS, RHS, mode);
    // PMOS_assigner(3, 1, 3, 4, 1, h, solution, LHS, RHS, mode);
    // NMOS_assigner(4, 4, 3, 0, 0, h, solution, LHS, RHS, mode);
    // PMOS_assigner(5, 1, 4, 5, 1, h, solution, LHS, RHS, mode);
    // NMOS_assigner(6, 5, 4, 0, 0, h, solution, LHS, RHS, mode);
}


// Voltage source stamp assigner
double Vs_assigner(int node_x, int node_y, double V_value, arma::mat &LHS, arma::mat &RHS){

    arma::vec value(1);
    value = V_value;
    // Extending the branch at the LHS matrix
    LHS = branch_ext(LHS, node_x, node_y);

    // Assigning the value at RHS
    RHS = arma::join_cols(RHS, value);

    // location of voltage value for transient simulation
    double V_locate1 = RHS.n_rows;

    return V_locate1;
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
    R_assigner(node_x,node_y,cond(x),LHS,RHS);
    // Matrix stamp for a capacitor on RHS
    Is_assigner(node_x,node_y,x1,LHS,RHS);
}

// Voltage pulse assigner
double V_pulse(double V1, double V2, double &t1,double td, double tr, double tf, double tpw, double tper, double h){
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
    t1 = t1 + h;
    return V_t1;
}

// Diode stamp assigner
void Diode_assigner(int node_x, int node_y, double Is, double VT, double cd, double h, arma::mat &LHS, arma::mat &RHS, arma::mat solution, int mode){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;

    double x = 0;
    double x1 = 0;
    double x2 = 0;
    double val_nodex = 0;
    double val_nodey = 0;

    

    if(mode == 0){
        cd = 0;
    }

    if((node_x == 0 && node_y == 0)){
        x = 0;
        x1 = 0;
    }
    else if(node_x == 0){
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
    
    // Matrix stamp for a diode on RHS
    Is_assigner(node_x,node_y,x1,LHS,RHS);
    // Matrix stamp for a diode on LHS
    R_assigner(node_x,node_y,cond(x),LHS,RHS);
}

// Newton Raphson system solver for non-linear and dynamic elements
arma::mat NewtonRaphson_system(arma::mat const init_LHS, arma::mat const init_RHS, arma::mat LHS, arma::mat RHS, arma::mat &solution, double h, int mode){
    int col_size = LHS.n_cols;
    int row_size = LHS.n_rows;
    double eps_val = 1e-8;
    arma::mat error = arma::zeros(1,1);
    double error_val = 9e9;
    error.row(0) = error_val;
    int iteration_counter = 0;
    arma::mat delta = arma::zeros(row_size,1);
    while((error(0,0) > eps_val) && (iteration_counter < 50)){ // iteration counter can be changed depending on the non-linearity of the circuit
        auto matrices = DynamicNonLinear(LHS,RHS,solution,h,mode);
        delta = arma::solve(matrices.first,(matrices.first*solution) - matrices.second);
        error.row(0) = arma::max(arma::abs(delta));
        solution -= delta;
        iteration_counter += 1;
    }
    return solution;
}

// Updates the value of RHS
arma::mat RHS_update(std::vector<double> RHS_locate, arma::mat RHS, std::vector<double> &val){
    for(int i = 0; i < RHS_locate.size(); i++){
        RHS.row(RHS_locate[i]-1).col(0) = val[i];
    }

    return RHS;
}

#endif