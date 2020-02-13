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

/** @file light_zxy_long_double.c
@brief memory functions for light zxy arrays
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


void malloc_light_zxy_long_double(struct dim_light *dim, long double * (***var))
{
	int x=0;
	int y=0;
	int z=0;

	*var = (long double ***) malloc(dim->zlen * sizeof(long double **));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (long double **) malloc(dim->xlen * sizeof(long double*));
		for (x = 0; x < dim->xlen; x++)
		{
			(*var)[z][x] = (long double *) malloc(dim->ylen * sizeof(long double));
			memset((*var)[z][x], 0, dim->ylen * sizeof(long double ));
		}
	}

}



void free_light_zxy_long_double(struct dim_light *dim, long double * (***in_var))
{
	int x=0;
	int y=0;
	int z=0;
	long double ***var=*in_var;

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

void flip_light_zxy_long_double_y(struct simulation *sim, struct dim_light *dim,long double *** data)
{
	int x=0;
	int y=0;
	int z=0;


	long double ***temp;

	malloc_light_zxy_long_double(dim, &temp);


	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				temp[z][x][y]=data[z][x][y];
			}
		}
	}


	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				data[z][x][dim->ylen-y-1]=temp[z][x][y];
			}
		}
	}


	free_light_zxy_long_double(dim, &temp);
}

void memset_light_zxy_long_double(struct dim_light *dim, long double ***data,int val)
{
	int x=0;
	int z=0;

	for (z = 0; z < dim->zlen; z++)
	{
		for (x = 0; x < dim->xlen; x++)
		{
			memset(data[z][x], val, dim->ylen * sizeof(long double ));
		}
	}

}

void div_light_zxy_long_double(struct dim_light *dim, long double ***data,long double val)
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
				data[z][x][y]/=val;
			}
		}
	}

}

//This shoudl be 3D interpolation but we are assuming the meshes are aligned.
long double interpolate_light_zxy_long_double(struct dim_light *dim, long double ***data,int z, int x, long double y_in)
{
	int y=0;
	long double x0=0.0;
	long double x1=0.0;
	long double y0=0.0;
	long double y1=0.0;

	long double ret;

	if (y_in<dim->y[0])
	{
		return 0.0;
	}


	if (y_in>=dim->y[dim->ylen-1])
	{
		//printf("here %Le %Le\n",y_in,dim->y[dim->ylen-1]);
		y=dim->ylen-1;
		x0=dim->y[y-1];
		x1=dim->y[y];
		y0=data[z][x][y-1];
		y1=data[z][x][y];

	}else
	{
		y=search(dim->y,dim->ylen,y_in);
		//printf("%d\n",y);
		x0=dim->y[y];
		x1=dim->y[y+1];

		y0=data[z][x][y];
		y1=data[z][x][y+1];
	}
	ret=y0+((y1-y0)/(x1-x0))*(y_in-x0);

return ret;

}

