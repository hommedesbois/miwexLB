/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

using namespace std;
#include "ludwig.h"

#ifndef __cplusplus
#error C++ compiler is required
#endif


extern "C" void DumpMacro(double **macro, int index){

    int off = current_slot * XMAXP * YMAXP;

	stringstream output_filename;
	output_filename << "vtk/FluidSimulationDomain_it_" << index << ".vtk";
	ofstream output_file;


	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS "<< NX << " " << NY << " " << 1 << endl;

	output_file << "X_COORDINATES " << NX << " float\n";
	for(int i=0; i<NX; i++)
		output_file << (double) i*DX << " ";
	output_file << endl;

	output_file << "Y_COORDINATES " << NY << " float\n";
	for(int j=0; j<NY; j++)
		output_file << (double) j*DX << " ";
	output_file << endl;

    output_file << "Z_COORDINATES " << 1 << " float\n";
    output_file << 0.01 << endl;


	output_file << "CELL_DATA " << XMAX * YMAX << endl;
	output_file << "SCALARS rho_prime float 1 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=YMAX; Y++)
			for(int X=1; X<=XMAX; X++){

				int idx = off + IDM(X,Y);
				output_file << (macro[idx][RHO] - rho0) << endl;
			}
    output_file << "SCALARS velocity float 3 \n";
	output_file << "LOOKUP_TABLE default\n";

		for(int Y=1; Y<=YMAX; Y++)
			for(int X=1; X<= XMAX; X++){

				int idx = off + IDM(X,Y);
				double ux = (macro[idx][RHOUX]/macro[idx][RHO])*Csound*sqrt(3);
				double uy = (macro[idx][RHOUY]/macro[idx][RHO])*Csound*sqrt(3);
				double uz = 0.0;
				output_file << ux << " " << uy << " " << uz << endl;
			}

	output_file.close();
}
