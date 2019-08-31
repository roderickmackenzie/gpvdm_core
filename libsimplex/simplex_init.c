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

/** @file simplex_init.c
@brief Init the simplex code
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <simplex.h>
#include <lock.h>
#include <const.h>

void multimin_init(struct multimin *data)
{
	data->ittr=0;
	data->n_max=100;
	data->ndim=0;
	data->nsimplex=data->ndim+1;
	data->stop_error=1e-9;
	data->p=NULL;
	data->center=NULL;
	data->y=NULL;
	data->fn=NULL;
	data->x=NULL;
	data->s=NULL;
	data->error=0.0;
	data->error_delta=1.0;
	data->error_last=1.0;
	if (lock_is_trial()==TRUE)
	{
		exit(0);
	}

}

void multimin_malloc(struct multimin *data)
{
	data->x=(double*)malloc(sizeof(double)*data->ndim);
	data->s=(double*)malloc(sizeof(double)*data->ndim);
	data->center=(double*)malloc(sizeof(double)*data->ndim);

	data->nsimplex=data->ndim+1;

	int s=0;
	int d=0;

	data->p=malloc(sizeof(double)*data->nsimplex);
	
	for (s=0;s<data->nsimplex;s++)
	{
		data->p[s]=malloc(sizeof(double)*data->ndim);
	}

	data->y=malloc(sizeof(double)*data->nsimplex);

	data->ptry=(double*)malloc(sizeof(double)*data->ndim);

	for (d=0;d<data->ndim;d++)
	{
		data->s[d]=0.0;
	}
}

void multimin_init_simplex(struct multimin *data)
{
	int d;
	int s;
	double delta=0.0;
	if (data->s[0]==0)
	{
		printf("Error: data->s[0]==0\n");
		return;
	}
	//data->p=malloc(sizeof(double)*data->nsimplex);
	
	for (s=0;s<data->nsimplex;s++)
	{

		for (d=0;d<data->ndim;d++)
		{
			delta=0.0;
			if (s==d)
			{
				delta=data->x[d]*data->s[d];
			}

			data->p[s][d]=data->x[d]+delta;
		}
	}


	for (s=0;s<data->nsimplex;s++)
	{

			data->y[s]=data->fn(data->p[s],data->ndim);
	}

	for (d=0;d<data->ndim;d++)
	{
		data->center[d]=-1;
	}

	//multimin_dump(data);
}

void multimin_free(struct multimin *data)
{
	int s;

	free(data->x);
	free(data->s);

	free(data->center);

	for (s=0;s<data->nsimplex;s++)
	{
		free(data->p[s]);
	}
	free(data->p);


	free(data->y);
	free(data->ptry);
}
