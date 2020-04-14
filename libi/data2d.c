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




/** @file data2d.c
	@brief Simple functions to read 2d scientific data from text files and perform simple maths on the data.
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sim_struct.h>

#include "data2d.h"
#include "util.h"
#include "cal_path.h"
#include "gpvdm_const.h"
#include <log.h>


void data2d_init(struct data2d *in, int x_len, int y_len)
{
	int i=0;
	in->y_len=y_len;
	in->x_len=x_len;
	in->y_mesh=NULL;
	in->x_mesh=NULL;

	in->data=(long double**)malloc(sizeof(long double*)*in->x_len);

	for (i=0;i<in->x_len;i++)
	{
		in->data[i]=(long double*)malloc(sizeof(long double)*in->y_len);
	}
}

void data2d_set_value(struct data2d *in,long double value)
{
	int x;
	int y;

	for (x=0;x<in->x_len;x++)
	{

		for (y=0;y<in->y_len;y++)
		{
			in->data[x][y]=value;
		}

	}
}

void data2d_free(struct data2d *in)
{
	int x;

	for (x=0;x<in->x_len;x++)
	{
		free(in->data[x]);
	}

	free(in->data);

	if (in->y_mesh!=NULL)
	{
		free(in->y_mesh);
	}


	if (in->x_mesh!=NULL)
	{
		free(in->x_mesh);
	}
}

void data2d_init_y_mesh(struct data2d *in,long double start, long double stop)
{
	int i;
	long double delta=(stop-start)/((long double)in->y_len);
	long double pos=start;
	in->y_mesh=(long double*)malloc(sizeof(long double)*in->y_len);

	for (i=0;i<in->y_len;i++)
	{
		pos+=delta;
		in->y_mesh[i]=pos;
	}

}


