/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include <cmath>
#include "pbe.h"
using namespace std;

void PBE::set_param()
{
  tend = 100.0;
  n_tstep =20;                                        // end time step
  v_min=1e-5;                                     // minimum size /volume of grid 
  v_max=100.0;                                      // maximum size / volume of grid
  r=1.5;                                          // geometric grid ratio
  m_v= ceil(log10(v_max/v_min)/log10(r)) +1 ;   // number of bin boundaries  
  m = m_v -1;                                    // number of pivots 
    
}

