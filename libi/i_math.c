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




/** @file i_mem.c
	@brief Memory management for i.c
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


static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));


/**Norm istruct to max value
@param in input istruct
*/
long double inter_norm(struct istruct* in,long double mul)
{
int i;
long double max=in->data[0];
//if (in->len>0) max=in->data[0];
for (i=0;i<in->len;i++)
{
if (in->data[i]>max) max=in->data[i];
}

for (i=0;i<in->len;i++)
{
in->data[i]*=mul/max;
}

return max;
}

/**Add together istruct structures
@param out output structure
@param in structure to add
*/
void inter_add(struct istruct* out,struct istruct* in)
{
int i;
	for (i=0;i<out->len;i++)
	{
		out->data[i]+=in->data[i];
	}

}

/**Get the average value of the data in a 1D istruct
@param in input istruct
*/
long double inter_avg(struct istruct* in)
{
int i;
long double sum=0.0;

for (i=0;i<in->len;i++)
{
	sum+=in->data[i];
}
return sum/((long double)(in->len));
}

/**Multiply the data in a 1D istruct by a number
@param in input istruct
*/
void inter_mul(struct istruct* in,long double mul)
{
int i;
//long double sum=0.0;

for (i=0;i<in->len;i++)
{
	in->data[i]*=mul;
}
}

/**Get maximum value of an istruct
@param in input istruct
*/
long double inter_get_fabs_max(struct istruct* in)
{
int i;
long double max=fabs(in->data[0]);
//if (in->len>0) max=in->data[0];
for (i=0;i<in->len;i++)
{
if (fabs(in->data[i])>max) max=fabs(in->data[i]);
}

return max;
}

/**Make all the data positive
@param in the structure holding the data
*/
void inter_mod(struct istruct* in)
{
int i;
for  (i=0;i<in->len;i++)
{
if (in->data[i]<0.0) in->data[i]*= -1.0;
}

}

/**Raise the data in an istruct by a power
@param p power to raise the data by
*/
void inter_pow(struct istruct* in,long double p)
{
int i;
for  (i=0;i<in->len;i++)
{
in->data[i]=pow(in->data[i],p);
}

}

/**Subtract two arrays they must be of the same length/x-asis
@param in opperand one, then result
@param in opperand two

*/
void inter_sub(struct simulation *sim,struct istruct* one,struct istruct* two)
{
if (one->len!=two->len)
{
	printf_log(sim,"The arrays are not the same length\n");
	exit(0);
}

int i;
for  (i=0;i<one->len;i++)
{
	if (one->x[i]!=two->x[i])
	{
		printf_log(sim,"The arrays do not have the same x axis\n");
		exit(0);
	}
	one->data[i]-=two->data[i];
}

}
