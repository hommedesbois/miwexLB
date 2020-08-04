/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <stdio.h>
#include "ludwig.h"

void BGKorDRTCollision(double *f, double **macro, int co){
int i, j, l, idf, idm, off;
double rho, rho_u, rho_v, u, v, feq, f_temp;
double Pxx_eq, Pxy_eq, Pyy_eq;
double Pxx_neq, Pxy_neq, Pyy_neq;

off = current_slot * xMaxp * yMaxp;

for(i=1; i<=xMax; i++)
    for(j=1; j<=yMax; j++){

        idm = off + IDM(i,j);

        rho = macro[idm][RHO];
        rho_u = macro[idm][RHOUX];
        rho_v = macro[idm][RHOUY];

        u = rho_u/rho;
        v = rho_v/rho;


        Pxx_eq = rho * (1./3. + (u * u));
        Pxy_eq = rho * u * v;
        Pyy_eq = rho * (1./3. + (v * v));

        Pxx_neq     = macro[idm][PXX]-Pxx_eq;
        Pxy_neq     = macro[idm][PXY]-Pxy_eq;
        Pyy_neq     = macro[idm][PYY]-Pyy_eq;

        for(l=0; l<NPOP; l++){

            idf = off*NPOP + IDF(i,j,l);

            feq = w[l]*rho*(1.0 -
            1.5*(u*u + v*v) +
            3.0*(ex[l]*u + ey[l]*v) +
            4.5*(ex[l]*u + ey[l]*v) *    (ex[l]*u + ey[l]*v) +
            0.5*(ex[l]*u + ey[l]*v) * (9*(ex[l]*u + ey[l]*v) * (ex[l]*u + ey[l]*v) - 9*(u*u+v*v)));

            switch(co){
            case 0:
                f[idf] = (1-omega)*f[idf] + omega * feq;
            break;
            case 1:
                f_temp = f[idf] - omega_n * (f[idf]-feq);
                f[idf] = f_temp + CoefTauNTau * w[l]*((ex[l]*ex[l]-Cs2)*Pxx_neq + 2*ex[l]*ey[l]*Pxy_neq + (ey[l]*ey[l]-Cs2)*Pyy_neq);
            break;
            }
        }
    }
}


void HRRCollision(double *cell, double **macro, double **grad){
    int i, j, l, idf, idm, off;
    double rho, rho_u, rho_v, u, v, feq, f_1;
    double Pxx_eq, Pxy_eq, Pyy_eq;
    double Pxx_neq, Pxy_neq, Pyy_neq;
    double Sxx, Syy, Sxy;
    double axx, ayy, axy;


    off = current_slot * xMaxp * yMaxp;


    for(i=1; i<=xMax; i++)
        for(j=1; j<=yMax; j++){


          idm = off + IDM(i,j);

          rho = macro[idm][0];
          rho_u = macro[idm][1];
          rho_v = macro[idm][2];

          u = rho_u/rho;
          v = rho_v/rho;

          //Sxx = grad[id_m][UXX] - 0.5 * (grad[id_m][UXX] + grad[id_m][UYY]);
          //Syy = grad[id_m][UYY] - 0.5 * (grad[id_m][UXX] + grad[id_m][UYY]);
          Sxx = grad[idm][UXX] ;
          Syy = grad[idm][UYY] ;
          Sxy = 0.5 * (grad[idm][UXY] + grad[idm][UYX]);

          Pxx_eq = rho * (1./3. + (u * u));
          Pxy_eq = rho * u * v;
          Pyy_eq = rho * (1./3. + (v * v));

          Pxx_neq     = macro[idm][3]-Pxx_eq;
          Pxy_neq     = macro[idm][4]-Pxy_eq;
          Pyy_neq     = macro[idm][6]-Pyy_eq;

          axx = Pxx_neq * sigma + (1.0-sigma)*(-2.0*tau*rho*Sxx*(1./3.)); 
          ayy = Pyy_neq * sigma + (1.0-sigma)*(-2.0*tau*rho*Syy*(1./3.)); 
          axy = Pxy_neq * sigma + (1.0-sigma)*(-2.0*tau*rho*Sxy*(1./3.)); 


          for(l=0; l<NPOP; l++){

            idf = off*NPOP + IDF(i,j,l);

            feq = w[l]*rho*(1.0 -
            1.5*(u*u + v*v) +
            3.0*(ex[l]*u + ey[l]*v) +
            4.5*(ex[l]*u + ey[l]*v) *    (ex[l]*u + ey[l]*v) +
            0.5*(ex[l]*u + ey[l]*v) * (9*(ex[l]*u + ey[l]*v) * (ex[l]*u + ey[l]*v) - 9*(u*u+v*v)));

            f_1 = w[l]*0.5*invCs4*((ex[l]*ex[l]-Cs2)*axx+2.0*ex[l]*ey[l]*axy+(ey[l]*ey[l]-Cs2)*ayy); // second order regularization
           
            cell[idf] = (1-omega)*f_1+feq;
            }
        }
}