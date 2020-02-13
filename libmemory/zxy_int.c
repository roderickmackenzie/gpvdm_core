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

/** @file zxy_int.c
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



void malloc_3d_int(struct dimensions *dim, int * (***var))
{
	int x=0;
	int y=0;
	int z=0;


	*var = (int ***) malloc(dim->zlen * sizeof(int **));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (int **) malloc(dim->xlen * sizeof(int*));
		for (x = 0; x < dim->xlen; x++)
		{
			(*var)[z][x] = (int *) malloc(dim->ylen * sizeof(int));
			memset((*var)[z][x], 0, dim->ylen * sizeof(int));
		}
	}

}

void free_3d_int(struct dimensions *dim, int ***var)
{
	int x=0;
	int y=0;
	int z=0;

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

}

void memory_flip_1d_int(int *var,int len)
{
	int x=0;
	int y=0;
	int z=0;
	int * data=malloc(sizeof(int)*len);
	for (y=0;y<len;y++)
	{
		data[y]=var[len-1-y];
	}

	for (y=0;y<len;y++)
	{
		var[y]=data[y];
	}

	free(data);

}

void dump_zxy_int(struct dimensions *dim, int ***var)
{
	int x=0;
	int y=0;
	int z=0;
	for (z = 0; z < dim->zlen; z++)
	{
		printf("z=%d:\n",z);

		for (y = 0; y < dim->ylen; y++)
		{

			for (x = 0; x < dim->xlen; x++)
			{
				printf("%d",var[z][x][y]);
			}
			printf("\n");
		}
	}
}
