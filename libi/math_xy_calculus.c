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




/** @file i.c
	@brief Simple functions to read in scientific data from text files and perform simple maths on the data.
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sim_struct.h>

#include "i.h"
#include "util.h"
#include "cal_path.h"
#include "gpvdm_const.h"
#include <log.h>
#include "inp.h"
#include <memory.h>


static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));

/**Integrate the data
@param in the structure to integrate
*/
long double inter_intergrate(struct math_xy* in)
{
int i;
long double tn=0.0;
long double tp=0.0;
//long double t=0.0;
long double dt=0.0;
//long double Eomega=0.0;
long double sum=0.0;
long double n;

	for (i=0;i<in->len;i++)
	{

		if (i==0)
		{
			tn=in->x[i];
		}else
		{
			tn=in->x[i-1];
		}

		if (i==in->len-1)
		{
			tp=in->x[i];
		}else
		{
			tp=in->x[i+1];
		}

		n=in->data[i];
		dt=fabs((tp-tn)/2.0);

		sum+=n*dt;


	}
return sum;
}

/**Integrate the data between limits
@param in the structure to integrate
@param from lower limit
@param from upper limit
*/
long double inter_intergrate_lim(struct math_xy* in,long double from, long double to)
{
int i;
long double tn=0.0;
long double tp=0.0;
//long double t=0.0;
long double dt=0.0;
//long double Eomega=0.0;
long double sum=0.0;
long double n;

	for (i=0;i<in->len;i++)
	{

		if (i==0)
		{
			tn=in->x[i];
		}else
		{
			tn=in->x[i-1];
		}

		if (i==in->len-1)
		{
			tp=in->x[i];
		}else
		{
			tp=in->x[i+1];
		}

		n=in->data[i];
		dt=fabs((tp-tn)/2.0);

		if (tn>from) sum+=n*dt;
		if (tn>to) break;

	}

return sum;
}



