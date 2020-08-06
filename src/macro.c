/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <stdio.h>
#include <math.h>
#include "ludwig.h"

/* Calculation of the macroscopic variables */

void ComputeMacroFromF(double *cell, double **macro){

int i, j, l, idf, idm;
double rho, rho_u, rho_v;
double Pxx, Pxy, Pyy;

int off = current_slot * XMAXP * YMAXP;

for(i=1; i<=XMAX; i++)
  for(j=1; j<=YMAX; j++){

      rho = 0.0, rho_u = 0.0, rho_v = 0.0;
      Pxx = 0.0, Pxy = 0.0, Pyy = 0.0;
      
      idm = off + IDM(i,j);
      
      for(l=0; l<NPOP; l++){

         idf = off * NPOP + IDF(i,j,l);

         rho         += cell[idf];
         rho_u       += cell[idf]*ex[l];
         rho_v       += cell[idf]*ey[l];
         Pxx         += cell[idf]*ex[l]*ex[l];
         Pxy         += cell[idf]*ex[l]*ey[l];
         Pyy         += cell[idf]*ey[l]*ey[l];
      }
   
        macro[idm][RHO] = rho;
        macro[idm][RHOUX] = rho_u;
        macro[idm][RHOUY] = rho_v;
        macro[idm][PXX] = Pxx;
        macro[idm][PXY] = Pxy;
        macro[idm][PYX] = Pxy;
        macro[idm][PYY] = Pyy;
    } // end i,j  loop
} // end function


void ComputeGradFromMacro(double **macro, double **grad){

int i, j, idnode, idxp, idxm, idyp, idym;

int off = current_slot * XMAXP * YMAXP;

for(i=1; i<=XMAX; i++)
  for(j=1; j<=YMAX; j++){

      idnode = off + IDM(i,j);
      idxp = off + IDM(i+1,j);
      idxm = off + IDM(i-1,j);
      idyp = off + IDM(i,j+1);
      idym = off + IDM(i,j-1);

      if(i==1)    
         idxm = off + IDM(XMAX,j);
      if(i==XMAX)
         idxp = off + IDM(1,j);
      if(j==1)
         idym = off + IDM(i, YMAX);
      if(j==YMAX)
         idyp = off + IDM(i,1);


      grad[idnode][UXX] = (macro[idxp][RHOUX]/macro[idxp][RHO] - macro[idxm][RHOUX]/macro[idxm][RHO]) * 0.5;
      grad[idnode][UYX] = (macro[idxp][RHOUY]/macro[idxp][RHO] - macro[idxm][RHOUY]/macro[idxm][RHO]) * 0.5;
      grad[idnode][UXY] = (macro[idyp][RHOUX]/macro[idyp][RHO] - macro[idym][RHOUX]/macro[idym][RHO]) * 0.5;
      grad[idnode][UYY] = (macro[idyp][RHOUY]/macro[idyp][RHO] - macro[idym][RHOUY]/macro[idym][RHO]) * 0.5;
      }
}
