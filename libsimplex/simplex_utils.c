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

/** @file simplex_utils.c
@brief Simplex utils
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <simplex.h>

void multimin_cal_center(struct multimin *data)
{
	#ifdef simplex_verbose
	printf("center\n");
	#endif
	//getchar();
	int s;
	int d;
	for (d=0;d<data->ndim;d++)
	{
		data->center[d]=0.0;
		for (s=0;s<data->nsimplex;s++)
		{
			if (s!=data->i_hi0)
			{
				data->center[d]+=data->p[s][d];
			}
		}

		data->center[d]/=(double)(data->nsimplex-1);
	}

}


void multimin_find_best(struct multimin *data)
{
	int s=0;
	int ss=0;
	
	int *ptrs=malloc(sizeof(int)*data->nsimplex);
	
	for (s=0;s<data->nsimplex;s++)
	{
		ptrs[s]=s;
	}

	int temp=0;
	int i=0;
	int swaped=1;
	while(swaped==1)
	{
		swaped=0;
		for (s=0;s<data->nsimplex-1;s++)
		{
			for (ss=0;ss<data->nsimplex;ss++)
			{
				if (data->y[ptrs[s]] < data->y[ptrs[s+1]])
				{
					temp=ptrs[s];
					ptrs[s]=ptrs[s+1];
					ptrs[s+1]=temp;
					swaped=1;
				}
			}
		}
	}
	//for (s=0;s<data->nsimplex;s++)
	//{
	//	printf("<%f> \n",data->y[ptrs[s]]);
	//}

	data->i_hi0=ptrs[0];
	data->i_hi1=ptrs[1];
	data->i_lo=ptrs[data->nsimplex-1];
	free(ptrs);
	
}

void sync(struct multimin *data,int s)
{
	int d=0;

	for (d=0;d<data->ndim;d++)
	{
		data->p[s][d] = data->ptry[d];
	}

	data->y[s]=data->ytry;
}

void simplex_copy_ans(struct multimin *data)
{
	int d;
	multimin_find_best(data);
	for (d=0;d<data->ndim;d++)
	{
		data->x[d]=data->p[data->i_lo][d];
	}

	data->error=data->y[data->i_lo];
}
/*
void swap(double *a,double *b)
{
	double temp=0;

	temp=*a;
	*a=*b;
	*b=temp;
}*/
