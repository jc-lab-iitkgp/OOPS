/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "nr3.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "odeint.h"
#include "quadrature.h"

class kernel
{
  public:
          double aggregation_kernel(double x, double y);        // aggregation kernel
          struct daughter_dist;                    // daughter distribution function
          struct Fun_1;
          struct Fun_2;
          double breakage_rate(double  x);                    // breakage rate equation
          double selection_func(double x);
};

struct kernel ::  daughter_dist{  // daughter distribution function
double v_prime;
daughter_dist(double v_prime) : v_prime(v_prime) {}
  Doub operator() (Doub v)        
    {  
        return (2/v_prime);
    }
};


struct kernel ::  Fun_1{          //Function used to find n_ij
double xk,xi,xj;
Fun_1(double xk,double xi,double xj) : xk(xk),xi(xi),xj(xj) {}
 
 Doub operator()(Doub x){
   return ((xj-x)*(2/xk)/(xj-xi));
 }
};

struct kernel ::  Fun_2{        //Function used to find n_ij
double xk,xi,xj;
Fun_2(double xk,double xi,double xj) : xk(xk),xi(xi),xj(xj) {}
 
 Doub operator()(Doub x){
   return ((x-xj)*(2/xk)/(xi-xj));
 }
};

                                 
            
