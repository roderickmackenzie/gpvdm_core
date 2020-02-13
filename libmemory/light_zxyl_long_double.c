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


void malloc_light_zxyl_long_double(struct dim_light *dim, long double * (****var))
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;

	*var = (long double ****) malloc(dim->zlen * sizeof(long double ***));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (long double ***) malloc(dim->xlen * sizeof(long double**));
		for (x = 0; x < dim->xlen; x++)
		{
			(*var)[z][x] = (long double **) malloc(dim->ylen * sizeof(long double*));
			for (y = 0; y < dim->ylen; y++)
			{
				(*var)[z][x][y] = (long double *) malloc(dim->llen * sizeof(long double));
				memset((*var)[z][x][y], 0, dim->llen * sizeof(long double ));

			}
		}
	}

}



void free_light_zxyl_long_double(struct dim_light *dim, long double * (****in_var))
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;

	long double ****var=*in_var;
	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zlen; z++)
	{

		for (x = 0; x < dim->xlen; x++)
		{
			for (y = 0; y < dim->ylen; y++)
			{
				free(var[z][x][y]);
			}
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

	*in_var=NULL;

}

void flip_light_zxyl_long_double_y(struct simulation *sim, struct dim_light *dim,long double **** data)
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;


	long double ****temp;

	malloc_light_zxyl_long_double(dim, &temp);

	for (l=0;l<dim->llen;l++)
	{
		for (z=0;z<dim->xlen;z++)
		{

			for (x=0;x<dim->xlen;x++)
			{

				for (y=0;y<dim->ylen;y++)
				{
					temp[z][x][y][l]=data[z][x][y][l];
				}

			}
		}
	}

	for (l=0;l<dim->llen;l++)
	{
		for (z=0;z<dim->zlen;z++)
		{

			for (x=0;x<dim->xlen;x++)
			{

				for (y=0;y<dim->ylen;y++)
				{
					data[z][x][dim->ylen-y-1][l]=temp[z][x][y][l];
				}
			}
		}
	}


	free_light_zxyl_long_double(dim, &temp);
}

void div_light_zxyl_long_double(struct dim_light *dim, long double ****data,long double val)
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;

	for (l = 0; l < dim->llen; l++)
	{
		for (z = 0; z < dim->zlen; z++)
		{
			for (x = 0; x < dim->xlen; x++)
			{
				for (y = 0; y < dim->ylen; y++)
				{
					data[z][x][y][l]/=val;
				}
			}
		}
	}

}

void memset_light_zxyl_long_double(struct dim_light *dim, long double ****data,int val)
{
	int x=0;
	int y=0;
	int z=0;

	for (z = 0; z < dim->zlen; z++)
	{
		for (x = 0; x < dim->xlen; x++)
		{
			for (y = 0; y < dim->ylen; y++)
			{
				memset(data[z][x][y], val, dim->llen * sizeof(long double ));
			}
		}
	}

}

void memset_light_zxyl_long_double_y(struct dim_light *dim, long double ****data,int z, int x, int l,long double val)
{
	int y=0;
	for (y = 0; y < dim->ylen; y++)
	{
		data[z][x][y][l]=val;
	}

}

