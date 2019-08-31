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

/** @file simplex_expand.c
@brief Simplex expand
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <simplex.h>

double expand(struct multimin *data,double mul)
{
	int d;

	multimin_cal_center(data);

	for (d=0;d<data->ndim;d++)
	{
		data->ptry[d]=data->center[d]+mul*(data->p[data->i_hi0][d]-data->center[d]);
	}

	data->ytry=(data->fn)(data->ptry,data->ndim);

	return data->ytry;
}

