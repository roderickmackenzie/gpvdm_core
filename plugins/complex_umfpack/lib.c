// 
// General-purpose Photovoltaic Device Model gpvdm.com- a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2012-2017 Roderick C. I. MacKenzie info at gpvdm dot com
// 
// https://www.gpvdm.com
// 
// 
// This program is free software; you can redistribute it and/or modify it
// under the terms and conditions of the GNU Lesser General Public License,
// version 2.1, as published by the Free Software Foundation.
// 
// This program is distributed in the hope it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
// 
// 

/** @file lib.c
	@brief Call the complex UMFPACK solver.
*/


#include <stdio.h>
#include <stdlib.h>
#include <lang.h>
#include <umfpack.h>
#include <util.h>

void complex_error_report(int status, const char *file, const char *func, int line)
{
	fprintf(stderr, "in %s: file %s, line %d: ", func, file, line);

	switch (status) {
		case UMFPACK_ERROR_out_of_memory:
			fprintf(stderr, _("out of memory!\n"));
			break;
		case UMFPACK_WARNING_singular_matrix:
			fprintf(stderr, _("matrix is singular!\n"));
			break;
		default:
			fprintf(stderr, "UMFPACK error code %d\n", status);
	}
}



int complex_umfpack_solver(struct simulation *sim,int c_col,int c_nz,int *c_Ti,int *c_Tj, long double *c_Tx,long double *c_Txz,long double *c_b,long double *c_bz)
{
	int i;
	void *Symbolic, *Numeric;
	int status;
	double *dtemp;
	int *itemp;

	if ((sim->c_last_col!=c_col)||(sim->c_last_nz!=c_nz))
	{
		dtemp = realloc(sim->c_x,c_col*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_x=dtemp;
		}

		dtemp = realloc(sim->c_xz,c_col*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_xz=dtemp;
		}

		dtemp = realloc(sim->c_b,c_col*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_b=dtemp;
		}

		dtemp = realloc(sim->c_bz,c_col*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_bz=dtemp;
		}

		itemp = realloc(sim->c_Ap,(c_col+1)*sizeof(int));
		if (itemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Ap=itemp;
		}

		itemp = realloc(sim->c_Ai,(c_nz)*sizeof(int));
		if (itemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Ai=itemp;
		}

		dtemp  = realloc(sim->c_Ax,(c_nz)*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Ax=dtemp;
		}

		dtemp  = realloc(sim->c_Az,(c_nz)*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Az=dtemp;
		}

		dtemp  = realloc(sim->c_Tx,(c_nz)*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Tx=dtemp;
		}


		dtemp  = realloc(sim->c_Txz,(c_nz)*sizeof(double));
		if (dtemp==NULL)
		{
			ewe(sim,"realloc failed\n");
		}else
		{
			sim->c_Txz=dtemp;
		}
		sim->c_last_col=c_col;
		sim->c_last_nz=c_nz;
	}

	for (i=0;i<c_col;i++)
	{
		sim->c_b[i]=(double)c_b[i];
		sim->c_bz[i]=(double)c_bz[i];

	}

	for (i=0;i<c_nz;i++)
	{
		sim->c_Tx[i]=(double)c_Tx[i];
		sim->c_Txz[i]=(double)c_Txz[i];
	}


	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];

	// get the default control parameters
	umfpack_zi_defaults (Control) ;

	//change the default print level for this demo
	//(otherwise, nothing will print)
	Control [UMFPACK_PRL] = 1 ;

	//print the license agreement
	//umfpack_zi_report_status (Control, UMFPACK_OK) ;
	Control [UMFPACK_PRL] = 0 ;

	// print the control parameters
	umfpack_zi_report_control (Control) ;

	status = umfpack_zi_triplet_to_col (c_col, c_col, c_nz, c_Ti, c_Tj, sim->c_Tx, sim->c_Txz, sim->c_Ap, sim->c_Ai, sim->c_Ax, sim->c_Az, NULL) ;


	if (status != UMFPACK_OK) {
		complex_error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

	// symbolic analysis
	//status = umfpack_di_symbolic(col, col, sim->complex_Ap, sim->complex_Ai, Ax, &Symbolic, NULL, NULL);
	status = umfpack_zi_symbolic(c_col, c_col, sim->c_Ap, sim->c_Ai, sim->c_Ax, sim->c_Az, &Symbolic, Control, Info) ;
	umfpack_zi_report_status (Control, status) ;

	if (status != UMFPACK_OK) {
		complex_error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

	// LU factorization
	//umfpack_di_numeric(sim->complex_Ap, sim->complex_Ai, sim->complex_Ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_zi_numeric (sim->c_Ap, sim->c_Ai, sim->c_Ax, sim->c_Az, Symbolic, &Numeric, Control, Info) ;

	if (status != UMFPACK_OK) {
		complex_error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}
	// solve system

	//umfpack_di_free_symbolic(&Symbolic);

		    // umfpack_di_solve(UMFPACK_A, sim->complex_Ap, sim->complex_Ai, sim->complex_Ax, x, b, Numeric, NULL, NULL);
	status = umfpack_zi_solve(UMFPACK_A, sim->c_Ap, sim->c_Ai, sim->c_Ax, sim->c_Az, sim->c_x, sim->c_xz, sim->c_b, sim->c_bz, Numeric, Control, Info) ;
	if (status != UMFPACK_OK) {
		complex_error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

		(void) umfpack_zi_report_vector (c_col, sim->c_x, sim->c_xz, Control) ;

	umfpack_zi_free_symbolic (&Symbolic) ;
	umfpack_di_free_numeric(&Numeric);

	for (i = 0; i < c_col; i++)
	{
		c_b[i]=sim->c_x[i];
		c_bz[i]=sim->c_xz[i];
	}

return 0;
}

void complex_umfpack_solver_free(struct simulation *sim)
{
	if (sim->c_x!=NULL)
	{
		free(sim->c_x);
		sim->c_x=NULL;	
	}

	if (sim->c_xz!=NULL)
	{
		free(sim->c_x);
		sim->c_x=NULL;	
	}

	if (sim->c_b!=NULL)
	{
		free(sim->c_b);
		sim->c_b=NULL;	
	}

	if (sim->c_bz!=NULL)
	{
		free(sim->c_b);
		sim->c_b=NULL;	
	}

	if (sim->c_Ap!=NULL)
	{
		free(sim->c_Ap);
		sim->c_Ap=NULL;
	}

	if (sim->c_Ai!=NULL)
	{
		free(sim->c_Ai);
		sim->c_Ai=NULL;
	}

	if (sim->c_Ax!=NULL)
	{
		free(sim->c_Ax);
		sim->c_Ax=NULL;
	}

	if (sim->c_Az!=NULL)
	{
		free(sim->c_Az);
		sim->c_Az=NULL;
	}

	if (sim->c_Tx!=NULL)
	{
		free(sim->c_Tx);
		sim->c_Tx=NULL;
	}

	if (sim->c_Txz!=NULL)
	{
		free(sim->c_Txz);
		sim->c_Txz=NULL;
	}
	sim->c_last_col=0;
	sim->c_last_nz=0;
}
