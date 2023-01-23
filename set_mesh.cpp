/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */ 

#include "pbe.h"
#include "nr3.h"

 
 void PBE::set_mesh()
 { 
   cout<<"Enter 0 for pure breakage,1 for pure aggregation,2 for aggregation+breakage by Fixed Pivot Technique and 3 for pure breakage,4 for pure aggregation,5 for aggregation+breakage by Finite Volume Technique"<<endl;
   wcin >>Process;

   double N= (log10(v_max/v_min)/log10(r)) +1;
   vi= new double[m_v];
   x= new double[m];
   w= new double [m];
   for( int i=0; i<m_v; i++)  vi[m_v-1-i]=v_min * pow(r,N-i-1); 

   for ( int i=0;i<m;i++)
     {
       x[i]=(vi[i+1] +vi[i])/2;
       w[i]=vi[i+1]-vi[i]; 
      
     }
    
 }
