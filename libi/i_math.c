//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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


/**Norm math_xy to max value
@param in input math_xy
*/
long double inter_norm(struct math_xy* in,long double mul)
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

/**Add together math_xy structures
@param out output structure
@param in structure to add
*/
void inter_add(struct math_xy* out,struct math_xy* in)
{
int i;
	for (i=0;i<out->len;i++)
	{
		out->data[i]+=in->data[i];
	}

}

/**Get the average value of the data in a 1D math_xy
@param in input math_xy
*/
long double inter_avg(struct math_xy* in)
{
int i;
long double sum=0.0;

for (i=0;i<in->len;i++)
{
	sum+=in->data[i];
}
return sum/((long double)(in->len));
}

/**Multiply the data in a 1D math_xy by a number
@param in input math_xy
*/
void inter_mul(struct math_xy* in,long double mul)
{
int i;
//long double sum=0.0;

for (i=0;i<in->len;i++)
{
	in->data[i]*=mul;
}
}

/**Get maximum value of an math_xy
@param in input math_xy
*/
long double inter_get_fabs_max(struct math_xy* in)
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
void inter_mod(struct math_xy* in)
{
int i;
for  (i=0;i<in->len;i++)
{
if (in->data[i]<0.0) in->data[i]*= -1.0;
}

}

/**Raise the data in an math_xy by a power
@param p power to raise the data by
*/
void inter_pow(struct math_xy* in,long double p)
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
void inter_sub(struct simulation *sim,struct math_xy* one,struct math_xy* two)
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
