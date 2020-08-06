/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <stdio.h>
#include <math.h>
#include "ludwig.h"

void InitializeFluid(double *cell, double **macro, int testcase){
    int i, j, l, id, idx;
    double rho, rho1,  ux,  uy;
    double Ppulse;
    double Es, Esc;

    double r0;
    double Rc = 0.1;

    double Xc = (double)(XMAX/2.)*DX;
    double Yc = (double)(YMAX/2.)*DX;

    double xpos, ypos; // position
    double r2;
    double Ueddy, Veddy;
    double duxdx, duxdy, duydx, duydy;

	double a0;
	double Length = DX * YMAX;

    double feq, fneq;

	for(i=0; i<XMAXP; i++)
		for(j=0; j<YMAXP; j++){

            id = IDM(i,j);
            // position of cell centres in 2D mesh with origin (0,0). First and last row of points lies outside the mesh --> required for periodicity condition
            xpos = (double) (i-1)*DX + (DX/2);
            ypos = (double) (j-1)*DX + (DX/2);

            switch(testcase){
                case 0:
                /*
                 ****** SHEAR LAYER *****
                 */
                r0 = 2.5*DX;
                rho = (P0/(Csound*Csound));
                a0 = 0.01*rho;

                U0 = 68.2; // Mach = 0.2
                //ux = U0*tanh((ypos+shift*DX-(Length/4.))/r0)*tanh((3.*Length/4.-(ypos+shift*DX))/r0);
                ux = U0*tanh((ypos-Length/4.)/r0)*tanh((3.*Length/4.-ypos)/r0);
		        uy = U0*a0*sin(2.*M_PI*(xpos/Length  + 0.25));

		        ux = ux / (sqrt(3.0)*Csound);
		        uy = uy / (sqrt(3.0)*Csound);
		        break;

                case 1:

		        /*
                 ****** PSEUDO ISENTROPIC VORTEX ****
                 */
                r2 = (xpos-Xc)*(xpos-Xc) + (ypos-Yc)*(ypos-Yc); // radius

                Es = 0.01; // vortex strength
                Esc = Es * Csound;
                
                rho1   = - ((rho0*Es*Es)/(2.0))*exp(1.0-(r2/(Rc*Rc)));
                rho     = rho0 + rho1 + 1./2.8 * (rho1 * rho1);

                U0 = 0.0; 
                V0 = 68.2; // Mach = 0.2
                
                Ueddy   = -(Esc/(Rc))*(ypos-Yc)*exp(0.5-(r2/(2.0*Rc*Rc)));
                ux      = (U0 + Ueddy) / (sqrt(3.0)*Csound);

                Veddy   = (Esc/(Rc))*(xpos-Xc)*exp(0.5-(r2/(2.0*Rc*Rc)));
                uy      = (V0 + Veddy) / (sqrt(3.0)*Csound);

                duxdx   =  Esc/Rc*(ypos-Yc)*exp(0.5-(r2/(2.0*Rc*Rc)))*(xpos-Xc)/(Rc*Rc);
                duxdx   *= DX / (sqrt(3.0)*Csound);

                duxdy   =  Esc/(Rc*Rc*Rc)*(ypos-Yc)*(ypos-Yc)*exp(0.5-(r2/(2.0*Rc*Rc))) - Esc/Rc * exp(0.5-(r2/(2.0*Rc*Rc)));
                duxdy   *= DX / (sqrt(3.0)*Csound);

                duydx   = -Esc/(Rc*Rc*Rc)*(xpos-Xc)*(xpos-Xc)*exp(0.5-(r2/(2.0*Rc*Rc))) + Esc/Rc * exp(0.5-(r2/(2.0*Rc*Rc)));
                duydx   *= DX / (sqrt(3.0)*Csound);

                duydy   = -Esc/Rc*(xpos-Xc)*exp(0.5-(r2/(2.0*Rc*Rc)))*(ypos-Yc)/(Rc*Rc);
                duydy   *= DX / (sqrt(3.0)*Csound);

                break;

                case 2:

                r0 = 0.002;
                Ppulse = 100.*exp(-((xpos-Xc)*(xpos-Xc)+(ypos-Yc)*(ypos-Yc))/r0);

                rho = (P0 + ((Pref+Ppulse) - Pref))/(Csound*Csound);

                U0 = 0.0; 
                V0 = 0.0;

                ux = U0 / (sqrt(3.0)*Csound); // in lattice units

                uy = V0 / (sqrt(3.0)*Csound);

                duxdx =  0.0;
                duxdx *= DX / (sqrt(3.0)*Csound);

                duxdy =  0.0;
                duxdy *= DX/ (sqrt(3.0)*Csound);

                duydx = 0.0;
                duydx *= DX/ (sqrt(3.0)*Csound);

                duydy = 0.0;
                duydy *= DX/ (sqrt(3.0)*Csound);

                break;
                }
            
            macro[id][0] = rho;
            macro[id][1] = rho * ux;
            macro[id][2] = rho * uy;

            for(l=0; l<NPOP; l++){
					idx = IDF(i,j,l);

                    feq = w[l]*rho*(1.0 -
						1.5*(ux*ux+uy*uy) +
						3.0*(ex[l]*ux+ey[l]*uy) +
						4.5*(ex[l]*ux+ey[l]*uy)*(ex[l]*ux+ey[l]*uy) +
                        0.5*(ex[l]*ux+ey[l]*uy)*(9*(ex[l]*ux+ey[l]*uy)*(ex[l]*ux+ey[l]*uy) - 9*(ux*ux+uy*uy)));


                    fneq = - w[l]*rho*tau*3.0*((ex[l]*ex[l]-Cs2)*duxdx + (ey[l]*ey[l]-Cs2)*duydy + ex[l]*ey[l]*(duxdy + duydx));

			cell[idx] = feq + fneq;
			// cell[idx] = feq;
        }
    }
}
