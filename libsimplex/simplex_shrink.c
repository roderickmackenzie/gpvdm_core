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

/** @file simplex_shrink.c
@brief Simplex shrink
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <simplex.h>

void multimin_shrink(struct multimin *data)
{
	printf("shrink\n");
int s=0;
int d=0;

	for (s=0;s<data->nsimplex;s++)
	{

		if (s != data->i_lo)
		{
			for (d=0;d<data->ndim;d++)
			{
				data->p[s][d]=0.5*(data->p[s][d]-data->p[data->i_lo][d]);
			}

			data->y[s]=(data->fn)(data->p[s],data->ndim);
		}
	}

	multimin_cal_center(data);
}


