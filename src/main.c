/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>
#include <time.h>

#include <sys/time.h>
#include <sched.h>

#include "ludwig.h" // lattice parameters

// constants model D3Q19
const double w[NPOP] = {4./9., // center
		1./9., 1./9., 1./9., 1./9., // face
		1./36., 1./36., 1./36., 1./36.}; //edge


/* discrete velocity
 *
 * index :            0   1   2   3   4   5   6   7   8
                      |   |   |   |   |   |   |   |   | */
const int ex[NPOP] = {0,  0,  0,  1, -1,  1, -1,  1, -1};
const int ey[NPOP] = {0,  1, -1,  0,  0,  1,  1, -1, -1};
/*                    |   |   |   |   |   |   |   |   | */
const int finv[NPOP]={0,  2,  1,  4,  3,  8,  7,  6,  5};


const double Csound = 343.2; // sound speed
// shared variables
double Mach; // mach number
double U0;
double V0;


double rho0 = 1.0; // density reference
double Pref; // pressure reference
double P0; // pressure parameter of solution -> p_solver

double viscosity = 1.5e-05; // viscosity

double DX;

// relaxation parameters
double omega, tau; 
double tau_n, omega_n;
double CoefTauNTau;

double Cs2 = 1./3.;
double invCs4 = 9.0;

int current_slot = 0, other_slot = 1; // two-field configuration

typedef struct idname{
    int index;
    char *name;
    } idname;

/** Initialization
  */
void InitializeFluid(double *, double **, int); 
/**  Propagate
  */
void Streaming(double *);
// NS
/**  Collision
  */
void HRRCollision(double *, double **, double **);
void BGKorDRTCollision(double *, double **, int);
/**  Macroscopic variables
  */
void ComputeMacroFromF(double *, double **);
void ComputeGradFromMacro(double **, double **);
/** Periodic BC
  */
void ComputeFcol_BC(double *);
/** Check conservation of mass
  */

void ComputeMass(double *);
/**  Write VTK output
  */
void DumpMacro(double **, int);


int main (int argc, const char * argv[]){
  int i;
  struct idname testcase, collision_operator;

  Pref = (Csound*Csound)/1.4;
  P0 = Csound * Csound;

  DX = 0.01;
  
  // non-dimensional relaxation time
  tau = 0.5 + sqrt(3)*viscosity/(Csound * DX);
  omega = 1. / tau;
  
  tau_n = 0.55;
  omega_n = 1./tau_n;

  CoefTauNTau = 0.5 * invCs4 * (tau-tau_n)/(tau*tau_n);

  // Allocation of space for distribution functions and macroscopic variables 
  
  double *f;
  double **macro, **grad;

  f = (double *) malloc (NPOP_TOTAL * sizeof(double));

  macro = (double **) malloc (2*XMAXP*YMAXP*sizeof(double *));
  for(i=0; i<(2 * XMAXP * YMAXP); i++){
        macro[i] = (double *)malloc(sizeof(double)*10);
        }

  grad = (double **) malloc (2*XMAXP*YMAXP*sizeof(double *));
  for(i=0; i<(2 * XMAXP * YMAXP); i++){
        grad[i] = (double *)malloc(sizeof(double)*4);
        }

/** Selcet collision model 
    *  0, BGK  
    *  1, DRT (Reference: https://doi.org/10.1103/PhysRevE.99.063305)
    *  2, HRR (Reference: https://doi.org/10.1080/14685248.2018.1540879)
    */    
  collision_operator.index=0;
  collision_operator.name="BGK";

/** Select testcase 
    *  0, shearlayer --> for doubly periodic shear layer 
    *  1, vortex --> for pseudo-isentropic vortex 
    *  2, pulse --> for acoutic pulse
    */ 
  testcase.index = 1;
  testcase.name = "vortex";

  fprintf(stdout, "#################################################################################\n");
  fprintf(stdout, "                                                                                #\n");
  fprintf(stdout, "       Simulation of the %s test case using the %s collision operator \n", testcase.name, collision_operator.name);
  fprintf(stdout, "                                                                                #\n");
  fprintf(stdout, "#################################################################################\n");
/** Simulation parameters
	  */
  fprintf(stdout, "# SIMULATION PARAMETERS #\n");
  fprintf(stdout, "Cs = %3.2f m/s \n", Csound);
  fprintf(stdout, "tau = %g \n", tau);
	fprintf(stdout, "viscosity = %g \n", viscosity);
  fprintf(stdout, "Dx = %g m \n", DX);
  fprintf(stdout, "Dt = %g s \n", (DX/(sqrt(3)*Csound)));
  
  
  fprintf(stdout, "#########################\n");


  // Initialize simulation
  InitializeFluid(f, macro, testcase.index); 

	int s;
	int index = 0;

	system("mkdir -p vtk");
  system("rm vtk/*.*");

  fprintf(stdout, "# ITERATIONS #\n");

	for(s=0; s<NSTEPS; s++){

/**  PROPAGATION
  */
    Streaming(f); // LBM

    other_slot = current_slot;
    current_slot = 1 - current_slot;

/**  COLLISION
  */
    ComputeMacroFromF(f, macro);
    switch(collision_operator.index){
      case 0:
        BGKorDRTCollision(f, macro, collision_operator.index);
      break;
      case 1:
        BGKorDRTCollision(f, macro, collision_operator.index);
      break;
      case 2: 
        ComputeGradFromMacro(macro, grad);
        HRRCollision(f, macro, grad);
      break;
    }
    
    ComputeFcol_BC(f);

    if(!(s % PERIOD)){

      fprintf(stdout, "time-step : %d\n", s);
      double physical_time = (s * DX)/(sqrt(3)*Csound);
      fprintf(stdout, "simulated time: %f s \n", physical_time);
      
      DumpMacro(macro, index);

      index++;
      }

    }
    return 0;
}
