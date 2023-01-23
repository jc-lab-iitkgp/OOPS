/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "kernel.h"
#include "mesh.h"

double kernel :: aggregation_kernel(double x, double y)
    {
        return  x ;                          // kernel function 
    }


double  kernel :: breakage_rate(double  x)        // breakage rate equation
    {    
        return  x * x;                         // power law breakge 
    } 


double kernel ::  selection_func(double x)       // selection function 
    {  
       return x ;
    }



