/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"

void  PBE :: compute_weight()   
{
   bw_b = new double[m];
   bw_d = new double[m]; 
   aw_b = new double *[m];
   for(int i=0;i<m;i++) aw_b[i] = new double[m];
   aw_d = new double *[m];  
   for(int i=0;i<m;i++) aw_d[i] = new double[m];      
   double sum1, sum2, sum3;
               
   for(int i=1;i<m;i++)
   {
     sum1=0.0;
     sum2=0.0;
     sum3=0.0;
     for(int k =0;k<=i;k++)
     {
       sum1+=bij[k][i];
       sum2+=((x[i]-x[k])*(bij[k][i]));       
       sum3+=(x[k]*bij[k][i]);
     }
       bw_b[i] = x[i] * (sum1 -1)/sum2;      
       bw_d[i] = bw_b[i] * sum3/x[i];   
   }
  double n,d;
  for(int i=0;i<m;i++)
  {
   for(int l=0;l<I.size();l++)
   {
    if(i == I[l][0]){
       n=x[I[l][1]]+x[I[l][2]];  
       d=2*x[I[l][0]]-n; 
       aw_b[I[l][1]][I[l][2]]=n/d;          // the birth-weigth term            
    }
   }
   for(int j=0;j<m;j++)
   {
    if(x[i]+x[j] <= v_max){
       for(int l=0;l<m;l++){
           if((i == I[l][2] && j == I[l][1]) || (i == I[l][1] && j == I[l][2])){
              n=x[I[l][0]];                             
              d=2*x[I[l][0]]-(x[i]+x[j]);         
              aw_d[i][j]=n/d;                  // the death-weigth term
           }
       }
    }
    else aw_d[i][j]=0;
   }
  } 
   
  
}
