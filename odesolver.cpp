/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */

#include "pbe.h"
#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

struct PBE :: breakage_aggregation
{   
    wchar_t Process;
    double *Sk, *bw_b, *bw_d, *Gam;
    double **bij,**n,**dgnl, **a, **aw_b, **aw_d;
    double ***eta;
    vector<vector<int>> I;

 breakage_aggregation(wchar_t Process,double ***eta,double **n,double **dgnl,double **bij,double *Sk,double **aw_b,double **aw_d,double *bw_b,double *bw_d,double *Gam,double **a,vector<vector<int> >I) : Process(Process), eta(eta), n(n), dgnl(dgnl), bij(bij), Sk(Sk), aw_b(aw_b), aw_d(aw_d), bw_b(bw_b), bw_d(bw_d), Gam(Gam), a(a), I(I) {}
 void  operator() (const Doub t, VecDoub_I &N, VecDoub_O &dNdt )    
  { 
    int m=N.size();
    double *FVS_b,*FVS_a,*FP_a,*FP_b;
    FVS_a= new double[m];
    FVS_b= new double[m];
    FP_a= new double[m];
    FP_b= new double[m];
    
    for (int i=0;i<m; i++)        //for breakage
     {  
      double birth=0;
      for(int j=i;j< m;j++) birth += n[i][j] * Gam[j] * N[j] ; 
       
      FP_b[i]= birth - Gam[i] * N[i] ;     // complete discretised pbe   
     } 

   for (int i=0;i<m;i++)                // for aggregation
   {  
     double birth=0;
     for(int j=0;j<m;j++)
      { 
       for (int k = j; k <m; k++) birth += (1- 0.5 * dgnl[j][k]) * eta[i][j][k] * N[j]* N[k]* a[j][k];
      }
     double sum_N=0;
     for ( int k=0; k<m; k++) sum_N+= N[k] * a[i][k];    
    
     FP_a[i]= (-1)*N[i]* sum_N + birth ;    
   }
 
  double sum=0.0,birth,death;
  
   death = 0.0; 
  for (int i=0;i<m;i++)
   { 
       death = Sk[i] * bw_d[i] * N[i]; 
        birth = 0.0;           
       for(int k = i;k<m;k++)  birth += ( bw_b[k]  *Sk[k] * N[k] * bij[i][k]);
       FVS_b[i]=birth-death; 

    }
   birth = 0.0;
   death = 0.0;
    
   for(int i=0;i<m;i++)
   {
    for(int l=0;l<I.size();l++)
    {
     sum=0.0;
     if(i == I[l][0]) sum+=(a[i][l] * aw_b[I[l][1]][I[l][2]] * N[I[l][1]] * N[I[l][2]]);     
    }
    
    birth = 0.5 *  sum;
    sum = 0.0;
    for(int j=0;j<m;j++) sum+=(a[i][j] * aw_d[i][j] * N[i] * N[j]); 
    death = sum;
   
    FVS_a[i]=(birth-death);
                     
   }   
       switch(Process){
       case '0':
              for ( int i =0;i<m;i++)  dNdt[i]= FP_b[i];
       break;
       case '1':
              for ( int i =0;i<m;i++)  dNdt[i]= FP_a[i];
       break;
       case '2':
              for ( int i =0;i<m;i++)  dNdt[i]= FP_b[i] + FP_a[i];
       break;
       case '3':
              for ( int i =0;i<m;i++)  dNdt[i]= FVS_b[i];
       break;
       case '4':
              for ( int i =0;i<m;i++)  dNdt[i]= FVS_a[i];
       break;
       case '5':
              for ( int i =0;i<m;i++)  dNdt[i]= FVS_b[i] + FVS_a[i];
       break;
       }
    }   
} ;




void PBE :: odesolver()
{  
    VecDoub Nstart2(m);
    
    for ( int i=0;i<m;i++) Nstart2[i]=Nstart1[i]; 

    double atol=1e-5, rtol=atol, h1=0.01, hmin=0.0, tstart=0.0;     //integration parameter 
     // this is for breakage and aggregation  odes solver simulataneous h1-> first guess step size hmin-> minimum allowed step size
     Output out(n_tstep-1);
     breakage_aggregation c(Process,eta,n,dgnl,bij,Sk,aw_b,aw_d,bw_b,bw_d,Gam,a,I);
     Odeint<StepperDopr5<breakage_aggregation> > ode(Nstart2,tstart,tend,atol,rtol,h1,hmin,out,c); 
     ode.integrate();
    
    N = new double *[out.count];
    for (int i=0;i<out.count;i++) N[i]=new double[m];
     
   for (int i=0;i<out.count;i++)
        {
        for(int j=0;j<m;j++)   N[i][j] = out.ysave[j][i];  
        }  
   
 }
 
