/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */
#include "pbe.h"

void  PBE :: compute_bij()        // bij
{    
  bij= new double *[m];
  for (int i = 0; i < m; i++) 
  {
    bij[i] = new double [m];
  }
  
   for(int i = 0 ;i<m;i++)
    { 
      for(int k = i ;k<m;k++)
      {
         if (k == i)
         {
          daughter_dist func(x[k]);
          bij[i][k]=qtrap(func,vi[i],x[i]);
         }
         else
         {
          daughter_dist func(x[k]);
          bij[i][k]=qtrap(func,vi[i],vi[i+1]);
         }
      } 
    }                      
}
