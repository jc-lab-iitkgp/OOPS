/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"

void PBE:: kronecker()
{
  //  dgnl= new double [m*m];
  dgnl= new double *[m];
  for (int i = 0; i < m; i++) {
    dgnl[i] = new double[m];
  }

  for (int row = 0; row < m; row++)
      {
          for (int col = 0; col < m; col++)
          {
              // Checking if row is equal to column
              if (row == col)
                  {
                    // dgnl[row*m + col]= 1;
                    dgnl[row][col]=1;
            
                  }
              else
              {
                // dgnl[row*m + col]= 0;
                dgnl[row][col]=0;

               
              }
          }
          
      }
    
}
