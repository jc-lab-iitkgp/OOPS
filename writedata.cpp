/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */


#include "pbe.h"
using namespace std;

void PBE :: writedata()
{
 file1.open ("PSD.txt");
 file1<<setprecision(6)<<fixed;
 for (int i=0;i<n_tstep;i++)
  {
   for (int j=0;j<m;j++) file1<<N[i][j]<<"  ";
   file1<<endl;
  } 
 file1.close();

}
