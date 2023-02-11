/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */


#include "pbe.h"
using namespace std;

void PBE :: writedata()
{
 /* -------------------------------------- Archita ------------------------*/
 file1.open ("PSD.txt");
 file1<<setprecision(6)<<fixed;
 for (int i=0;i<n_tstep;i++)
  {
   for (int j=0;j<m;j++) file1<<N[i][j]<<"  ";
   file1<<endl;
  } 
 file1.close();
 /* ---------------------------------------------------------------------------
 ***********************  change Made by Manish Saini  **********************
 ---------------------------------------------------------------------------*/

 /* --------------------------------------------------------------------------- 
     *************** numerical Number density  distribution *************************
     ----------------------------------------------------------------------------*/
    ofstream hout_0("numerical_number_density_moment.txt");
     for (int i = 0; i < count; i++)
    {    if (i==85)
        {for (int  j = 0; j < m; j++)
        {    
            hout_0 << x[j]<<"   "<< N[j][i]/ w[j]<<endl;
        } }  
    }
    
    /* --------------------------------------------------------------------------- 
     *************** numerical zeroth moment distribution *************************
     ----------------------------------------------------------------------------*/
    ofstream hout("numerical_zeroth_moment.txt");
    double mu=0;
    for (int i = 0; i < count; i++)
    {
        double Ntot=0;
        for (int  j = 0; j < m; j++)
        {
            Ntot+= N[j][i];
        }
        if (i==0)
        {
             mu=Ntot;
        }

        hout<<t_Afsol[i] <<"  "<< Ntot/mu <<endl; 
            
    }
    hout.close();
    /* --------------------------------------------------------------------------- 
     *************** numerical first moment distribution *************************
     ----------------------------------------------------------------------------*/
    ofstream hout1("numerical_first_moment.txt");
    double mu1_t=0;
    for (int i = 0; i < count; i++)
    {
        double mu1=0;
        for (int  j = 0; j < m; j++)
        {
            mu1+= x[j] * N[j][i];
        }
        if (i==0) mu1_t=mu1;
        hout1<<t_Afsol[i] <<"  "<< mu1/mu1_t <<endl;   
    }
    hout1.close();
     /* --------------------------------------------------------------------------- 
     *************** numerical second moment distribution *************************
     ----------------------------------------------------------------------------*/
    ofstream hout2("numerical_second_moment.txt");
    for (int i = 0; i < count; i++)
    {
        double mu2=0;
        for (int  j = 0; j < m; j++)
        {
            mu2+= x[i]*x[i] * N[j][i];
        }
        hout2<<t_Afsol[i] <<"  "<< mu2 <<endl;  // y(t)= m*t +c m-> slope c-> intercept   
    }
    hout2.close();
    /* --------------------------------------------------------------------------- 
     *************** analytical zeroth moment distribution for pure breakage *****
     ----------------------------------------------------------------------------*/
    ofstream hout3("Pure_breakage_analytical_zeroth_moment.txt"); 
    for (int i = 0; i < count; i++)
    {
        double mu1 = 0;
        for ( int j=0;j< m;j++)
        {
            mu1+= x[j]* N[j][i];
        }
        hout3 <<t_Afsol[i]<<"  "<< (1 +mu1*t_Afsol[i])<<endl;   // y(t)= mt +c      
    }
    hout3.close();

    /* --------------------------------------------------------------------------- 
     *************** analytical zeroth moment distribution for pure breakage *****
     ----------------------------------------------------------------------------*/
    ofstream hout4("Pure_aggregation_analytical_zeroth_moment.txt");
    for (int i = 0; i < count; i++)
    {
        // for constant kernel 1 is multiply for value of the constant kernel
        hout4<< t_Afsol[i]<<"  "<< (mu/(1/mu + t_Afsol[i]/2))<<endl;  // mu is intial total number of particle 
    }
    hout4.close();

    /* --------------------------------------------------------------------------- 
     *************** Analytical  Number density  distribution *************************
     ----------------------------------------------------------------------------*/
    ofstream hout5("analytical_number_density_moment.txt");
     for (int i = 0; i < count; i++)
    {    
        if (i==85)
        {
        for (int  j = 0; j < m; j++)
        {    
            hout5 << x[j]<<"   "<< pow(1+t_Afsol[i], 2)* exp(- x[j]*(1+t_Afsol[i])) <<endl;
        } 
        }  
    }
   /*  ---------------------------------------------------------------------------
       *********************       end ***********************************
       ------------------------------------------------------------------------*/

}
