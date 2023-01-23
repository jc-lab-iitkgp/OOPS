/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"
#include <iostream>
void  PBE::compute_eta_ijk()
{
  eta=new double **[m];
  for (int i = 0; i < m; i++)
   {
        eta[i] = new double *[m];
        for (int j = 0; j < m; j++)
        {
            eta[i][j] = new double [m];
        }
    }

  for ( int j=0;j<m;j++)
      {
          for ( int k=j;k<m;k++)
          {   
              double v=0;
              v= x[j]+ x[k];
              for( int i=0; i<m-1;i++)
              { 
                  if(v>=x[i] && v<= x[i+1])
                  {
                     eta[i+1][j][k]= (v-x[i])/(x[i+1]-x[i]);
                     eta[i][j][k]= (x[i+1]-v)/(x[i+1]-x[i]);
                     break;

                  }
               
              }
          }
      }  
     
      
}

