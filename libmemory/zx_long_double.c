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

/** @file zx_long_double.c
@brief memory functions for 3D arrays
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include "memory.h"




void malloc_zx_gdouble(struct dimensions *dim, gdouble * (**var))
{
	int z=0;

	*var = (gdouble **) malloc(dim->zlen * sizeof(gdouble *));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (gdouble *) malloc(dim->xlen * sizeof(gdouble));
		memset((*var)[z], 0, dim->xlen * sizeof(gdouble));
	}

}

void mem_set_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in)
{
	int z=0;
	int x=0;


	for (z = 0; z < dim->zlen; z++)
	{
		for (x = 0; x < dim->xlen; x++)
		{
			data_out[z][x]=data_in[z][x];
		}
	}

}

void mem_add_zx_gdouble_from_zx_gdouble(struct dimensions *dim, gdouble **data_out, gdouble **data_in)
{
	int z=0;
	int x=0;


	for (z = 0; z < dim->zlen; z++)
	{
		for (x = 0; x < dim->xlen; x++)
		{
			data_out[z][x]+=data_in[z][x];
		}
	}

}

void free_zx_gdouble(struct dimensions *dim, gdouble * (**in_var))
{
	int z=0;
	long double **var=*in_var;

	if (var==NULL)
	{
		return;
	}

	//printf("%d %d %d\n",dim->zlen,dim->xlen,dim->ylen);
	for (z = 0; z < dim->zlen; z++)
	{
		free(var[z]);
	}

	//printf("free\n");
	free(var);

	*in_var=NULL;
}

void malloc_zx_int(struct dimensions *dim, int * (**var))
{
	int z=0;

	*var = (int **) malloc(dim->zlen * sizeof(int *));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (int *) malloc(dim->xlen * sizeof(int));
		memset((*var)[z], 0, dim->xlen * sizeof(int));
	}

}

void free_zx_int(struct dimensions *dim, int *(**in_var))
{
	int z=0;

	int **var=*in_var;

	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zlen; z++)
	{
		free(var[z]);
	}

	free(var);

	*in_var=NULL;
}

