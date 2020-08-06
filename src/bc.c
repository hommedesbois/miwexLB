/* T. Horstmann
* AT-TRA
* DLR
* Berlin
* Germany
*/

#include <stdio.h>
#include <math.h>
#include "ludwig.h"


void ComputeFcol_BC(double *f){
	int i,j,l;

	/*
	 *		Y
	 *		^
	 *		|
	 *		|
	 *		Z ----> X
	 */

	int i_in, i_out, j_in, j_out;
	int idx_in, idx_out, off;

	off = current_slot*XMAXP*YMAXP*NPOP;

	// faces
	i_out = 1;
	i_in = XMAX + 1;

	for(j=1; j<=YMAX; j++)
			for(l=0; l<NPOP; l++){
				idx_out =  off + IDF(i_out,j,l);
				idx_in = off + IDF(i_in,j,l);
				f[idx_in] = f[idx_out];
			}

	i_out = XMAX;
	i_in = 0;

	for(j=1; j<=YMAX; j++)
			for(l=0; l<NPOP; l++){
				idx_out =  off + IDF(i_out,j,l);
				idx_in = off + IDF(i_in,j,l);
				f[idx_in] = f[idx_out];
			}

	j_out = 1;
	j_in = YMAX + 1;

	for(i=1; i<=XMAX; i++)
    	for(l=0; l<NPOP; l++){
			idx_out =  off + IDF(i,j_out,l);
			idx_in = off + IDF(i,j_in,l);
			f[idx_in] = f[idx_out];
			}

	j_out = YMAX;
	j_in = 0;

	for(i=1; i<=XMAX; i++)
		for(l=0; l<NPOP; l++){
			idx_out =  off + IDF(i,j_out,l);
			idx_in = off + IDF(i,j_in,l);
			f[idx_in] = f[idx_out];
			}

	//edges
    i_out = 1;
    j_out = 1;
    i_in = XMAX + 1;
    j_in = YMAX + 1;
    for(l=0; l<NPOP; l++){
    	idx_out = off + IDF(i_out, j_out,l);
    	idx_in = off + IDF(i_in,j_in,l);
    	f[idx_in] = f[idx_out];
    }

    i_out = XMAX;
    j_out = YMAX;
    i_in = 0;
    j_in = 0;
    for(l=0; l<NPOP; l++){
    	idx_out = off + IDF(i_out, j_out,l);
    	idx_in = off + IDF(i_in,j_in,l);
    	f[idx_in] = f[idx_out];
    }

    i_out = XMAX;
    j_out = 1;
    i_in = 0;
    j_in = YMAX + 1;
    for(l=0; l<NPOP; l++){
    	idx_out = off + IDF(i_out,j_out,l);
    	idx_in = off + IDF(i_in,j_in,l);
    	f[idx_in] = f[idx_out];
    }

    i_out = 1;
    j_out = YMAX;
    i_in = XMAX + 1;
    j_in = 0;
    for(l=0; l<NPOP; l++){
    	idx_out = off + IDF(i_out,j_out,l);
    	idx_in = off + IDF(i_in,j_in,l);
    	f[idx_in] = f[idx_out];
    }
}