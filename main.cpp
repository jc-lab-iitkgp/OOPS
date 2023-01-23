/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"
using namespace std;


int main()
{   
    PBE pbe;
    pbe.set_param();              // setting the parameter
    pbe.set_mesh();               // creation of geometric  mesh
    pbe.initial_condition();      // setting the iniital condition for all nodes 
    pbe.kronecker();              // kronecker delta  matrics
    pbe.compute_eta_ijk();        // computing the eta_ijk matrix 
    pbe.kernel_a_ij();            // computing then a_ij matrix
    pbe.set_indices();            //Set formation for aggregation for FVS 
    pbe.compute_Sk();             //selection function
    pbe.compute_bij();            //integral value of daughter dist func 
    pbe.compute_weight();         //weight calculation
    pbe.compute_n_ij();           // computing the n_ij  matrix    
    pbe.odesolver();              // odesolver  using odeint.h
    pbe.writedata();    
    return 0;
}
