/*Code to solve 1D Population Balance Equation using Fixed Pivot Technique or Finite Volume Method .
Reference - [1]Jayanta Chakraborty , 'Engineering of Submicron Particles: Fundamental Concepts and Models' , 2019,John Wiley & Sons Ltd.

[2] Kumar et al.'An accurate and efficient discrete formulation of aggregation population balance equation', 2016, Kinetic & Related Models, vol-9, pp. 373-391
Developed by JC , Archita Karar and Manish Saini */


class mesh
{
    public:
          wchar_t Process;
          double atol;                 // absolute tolerence 
          double rtol;                 // relative tolerence
          double h1;                   // first guess step size 
          double hmin;                 // minimum allowed guess step size
          double tstart;               // start time step
          double tend;                 // end time step
          int n_tstep;              //no of time-steps
          double v_min;                // minimum volume/ size of grid
          double v_max;                // maximum volume / size of grid
          double r;                    // size ratio    
          int    m_v;                  // mumber of bin boundaries
          int    m;                    // number of volume nodes
          double *vi;                  // the bin boudaries  vector for 
          double *w;                   // bin width array
          double *x;                   // pivots array
          double *Nstart1;             // initial condition for all pivots
          double **N;
          double *Gam;                 // Gam for breakage rate 
          double ***eta;               // eta_ijk generation matrix for aggreagtion
          double **n;                  // n_ij generation matrix for breakage
          double **a;                  // kernel matrix
          double **dgnl;               // dgnl kronecker delta function
          vector<vector<int>> I;       //set of indices for aggregation
          double *bw_b,*bw_d;          //weights for breakage
          double **aw_b,**aw_d;          //weights for aggregation
          double **bij;                //integral of daughter distribution function
          double *Sk;                  //Selection function value
          double *t_Afsol;              // time step after solution
          double **Number_dist;     // number distribution after solution for each time step
          ofstream file1;
};
