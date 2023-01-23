/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"

void PBE :: set_indices()
{
  int n_tuple=0;
  
  I.resize(1,vector<int>(3));
  
  for(int i=0;i<m;i++)
  {
   for(int j=0;j<m;j++)
   {
    for(int k=0;k<m;k++)
    {
   
     if((x[j]+x[k]<=vi[i+1])  &&  (x[j]+x[k]>vi[i])){  
     
     I[n_tuple][0]=i;
     I[n_tuple][1]=j;
     I[n_tuple][2]=k;
     
     n_tuple++;
     I.resize(n_tuple+1,vector<int>(3));
     }
    }   
   }
  }
  I.resize(n_tuple,vector<int>(3));                             
}    
