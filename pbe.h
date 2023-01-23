/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "kernel.h"
#include "mesh.h"

class PBE: public kernel, public mesh
{ 
  public: 
          struct breakage_aggregation;
          void  set_param();           // setting the parameter
          void  set_mesh();            // creation of geometric  mesh 
          void  compute_eta_ijk();     // computing the eta_ijk matrix 
          void  kernel_a_ij();         // computing then a_ij matrix  
          void  initial_condition();   // setting the iniital condition for all nodes 
          void  compute_n_ij();        // computing the n_ij  matrix
          void  set_indices();
          void  compute_Sk();          //storing selection function value in Sk
          void  compute_bij();         //intral value of daughter dist function
          void  compute_weight();      //weight for breakage-aggregation
          void  kronecker();           // kronecker delta  matrics 
          void  odesolver();           // odesolver using odeint.h 
          void  writedata();
};

