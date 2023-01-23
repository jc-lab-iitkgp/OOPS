/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"

void PBE :: compute_n_ij()
{   
  n= new double *[m];
  for (int i = 0; i < m; i++) 
  {
    n[i] = new double [m];
  }
    for ( int i=0;i<m;i++)
    {
     if(i ==0)
      {
       for (int k=1;k<m;k++)
        { 
         Fun_1 func(x[k],x[i],x[i+1]);  
         n[i][k] = qtrap(func,x[i],x[i+1]) ; 
        }
       }
       else
       {
         Fun_2 func(x[i],x[i],x[i-1]);
         n[i][i] = qtrap(func,x[i-1],x[i]) ;
        for (int k=i+1;k<m;k++)
        {
         Fun_1 func1(x[i],x[i],x[i-1]);
         n[i][k] = qtrap(func,x[i],x[i+1]) ;

         Fun_2 func2(x[k],x[i],x[i-1]);
         n[i][k] += qtrap(func,x[i-1],x[i]) ;
        }
       }
    }
}


