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

/** @file solver.c
	@brief Complex UMFPACK solver interface.
*/


#include <solver_interface.h>
#include <dll_export.h>
#include <util.h>
#include "lib.h"
#include <log.h>


EXPORT void dll_complex_matrix_init(struct simulation *sim)
{
//printf("init\n");
}


EXPORT void dll_complex_matrix_solve(struct simulation *sim,int c_col,int c_nz,int *c_Ti,int *c_Tj, long double *c_Tx,long double *c_Txz,long double *c_b,long double *c_bz)
{
	complex_umfpack_solver(sim,c_col,c_nz,c_Ti,c_Tj, c_Tx,c_Txz,c_b,c_bz);
}

EXPORT void dll_complex_matrix_solver_free(struct simulation *sim)
{
	complex_umfpack_solver_free(sim);
}

