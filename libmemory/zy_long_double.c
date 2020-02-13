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

/** @file memory_basic.c
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


void malloc_zy_long_double(struct dimensions *dim, long double * (**var))
{
	int z=0;

	*var = (gdouble **) malloc(dim->zlen * sizeof(gdouble *));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (gdouble *) malloc(dim->ylen * sizeof(gdouble));
		memset((*var)[z], 0, dim->ylen * sizeof(gdouble));
	}

}

void free_zy_long_double(struct dimensions *dim, long double * (**in_var))
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

	free(var);

	*in_var=NULL;
}

void malloc_zy_int(struct dimensions *dim, int * (**var))
{
	int z=0;

	*var = (int **) malloc(dim->zlen * sizeof(int *));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (int *) malloc(dim->ylen * sizeof(int));
		memset((*var)[z], 0, dim->ylen * sizeof(int));
	}

}

void free_zy_int(struct dimensions *dim, int *(**in_var))
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

