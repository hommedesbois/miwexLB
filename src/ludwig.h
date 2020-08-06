/*
 *  ludwig.h
 *
 *  Created by  Tobias Horstmann
 *  Copyright 2020 DLR. All rights reserved.
 *
 */

#define NSTEPS 1001
#define PERIOD 100
#define SIGMA 0.98 // Coef for HRR scheme

enum {UXX=0, UXY, UYX, UYY};
enum {RHO=0, RHOUX, RHOUY, PXX, PXY, PYX, PYY};

#define max(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })

/*
 *
 *		Y
 *		^
 *		|
 *		|
 *		Z ----> X
 *
 *  periodic BC in all directions for root grid
 */

/*
 * DATA POINTS
 */

#define XMAX 100
#define YMAX 100

#define XMAXP (XMAX + 2) // padded dimensions for periodic BC and transition overlay
#define YMAXP (YMAX + 2)

/*
 * MESH
 */

#define NX (XMAX + 1)
#define NY (YMAX + 1)

/*
 * NUMBER OF DISTRIBUTIONS
 */ 

#define NPOP 9 // D2Q9 lattice
#define NPOP_TOTAL XMAXP * YMAXP * NPOP * 2

/*
 * INDEX FOR ARRAYS OF DISTRIBUTION AND MACROSCOPIC VARIABLES
 */  

#define IDF(i,j,l)  (NPOP*((j)+YMAXP*(i))+(l))
#define IDM(i,j)    (i) * YMAXP + (j)

// shared variables
extern const double w[NPOP];
extern const int ex[NPOP];
extern const int ey[NPOP];
extern const int finv[NPOP];

extern int current_slot, other_slot;


extern double rho0; //  density parameter
extern const double Csound; // sound speed
extern double U0, V0; //  velocity parameter
extern double Pref; // pressure reference
extern double P0; // pressure parameter
extern double viscosity; // viscosity
extern double Mach;

extern double Cs2;
extern double invCs4;

extern double Length;
extern double pos;
extern double DX;

extern double omega, tau, omega_n, tau_n;
extern double CoefTauNTau;
