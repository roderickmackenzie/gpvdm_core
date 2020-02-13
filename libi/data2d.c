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


