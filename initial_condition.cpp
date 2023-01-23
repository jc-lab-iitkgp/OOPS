/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"

struct initial_fun_integrator    
{
    double a,b;
    initial_fun_integrator( double a,double b) : a(a), b(b) {}
    double operator ()( double  x) 
    {   
        
         return  (1/0.01)*exp(-x/0.01);                 // function for initial condition on the nodes
    }
};


 void PBE :: initial_condition()
{   
    Nstart1 = new double[m];
    Gam = new double[m];
    for (int i=0;i<m;i++)
     {
      //initial_fun_integrator f(0,0);
      //Nstart1[i] = qtrap(f,vi[i], vi[i+1]);
      Nstart1[i]=0.0; 
      Gam[i]=breakage_rate(x[i]);    
     } 
     Nstart1[m-1]=1.0;   
}
