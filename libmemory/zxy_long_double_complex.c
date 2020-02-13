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


void malloc_zxy_long_double_complex(struct dimensions *dim, long double complex * (***var))
{
	int x=0;
	int y=0;
	int z=0;


	*var = (long double complex ***) malloc(dim->zlen * sizeof(long double complex **));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (long double complex **) malloc(dim->xlen * sizeof(long double complex*));
		for (x = 0; x < dim->xlen; x++)
		{
			(*var)[z][x] = (long double complex *) malloc(dim->ylen * sizeof(long double complex));
			memset((*var)[z][x], 0, dim->ylen * sizeof(long double complex));
		}
	}

}



void free_zxy_long_double_complex(struct dimensions *dim, long double complex * (***in_var))
{
	int x=0;
	int y=0;
	int z=0;

	long double complex ***var=*in_var;
	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zlen; z++)
	{

		for (x = 0; x < dim->xlen; x++)
		{
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

	*in_var=NULL;

}

