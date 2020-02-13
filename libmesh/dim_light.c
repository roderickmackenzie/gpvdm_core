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

/** @file dim.c
@brief Dimension object
*/

#include <device.h>
#include "mesh.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include "hard_limit.h"
#include <log.h>
#include <cal_path.h>
#include <lang.h>
#include <shape.h>

void dim_light_init_xyzl(struct dim_light *dim,char xyzl)
{
	if (xyzl=='x')
	{
		dim->x= NULL;
		dim->dx= -1.0;
		dim->xlen=-1;
	}else
	if (xyzl=='y')
	{
		dim->y= NULL;
		dim->dy= -1.0;
		dim->ylen=-1;
	}else
	if (xyzl=='z')
	{
		dim->z= NULL;
		dim->dz= -1.0;
		dim->zlen=-1;
	}else
	if (xyzl=='l')
	{
		dim->l= NULL;
		dim->dl= -1.0;
		dim->llen=-1;
	}

}

void dim_light_init(struct dim_light *dim)
{
	dim_light_init_xyzl(dim,'x');
	dim_light_init_xyzl(dim,'y');
	dim_light_init_xyzl(dim,'z');
	dim_light_init_xyzl(dim,'l');
}

void dim_light_free_xyzl(struct dim_light *dim,char xyzl)
{
	if (xyzl=='x')
	{
		if (dim->x!=NULL)
		{
			free(dim->x);
			dim_light_init_xyzl(dim,'x');
		}
	}else
	if (xyzl=='y')
	{
		if (dim->y!=NULL)
		{
			free(dim->y);
			dim_light_init_xyzl(dim,'y');
		}
	}else
	if (xyzl=='z')
	{
		if (dim->z!=NULL)
		{
			free(dim->z);
			dim_light_init_xyzl(dim,'z');
		}
	}else
	if (xyzl=='l')
	{
		if (dim->l!=NULL)
		{
			free(dim->l);
			dim_light_init_xyzl(dim,'l');
		}
	}

}

void dim_light_free(struct dim_light *dim)
{
	dim_light_free_xyzl(dim,'x');
	dim_light_free_xyzl(dim,'y');
	dim_light_free_xyzl(dim,'z');
	dim_light_free_xyzl(dim,'z');
	dim_light_init(dim);
}

void dim_light_malloc_xyzl(struct dim_light *dim,char xyzl)
{

	if (xyzl=='x')
	{
		dim->x = (long double *) malloc(dim->xlen * sizeof(long double));
		memset(dim->x, 0, dim->xlen * sizeof(long double));
	}else
	if (xyzl=='y')
	{
		dim->y = (long double *) malloc(dim->ylen * sizeof(long double));
		memset(dim->y, 0, dim->ylen * sizeof(long double));
	}else
	if (xyzl=='z')
	{
		dim->z = (long double *) malloc(dim->zlen * sizeof(long double));
		memset(dim->z, 0, dim->zlen * sizeof(long double));
	}else
	if (xyzl=='l')
	{
		dim->l = (long double *) malloc(dim->llen * sizeof(long double));
		memset(dim->l, 0, dim->llen * sizeof(long double));
	}


}



void dim_light_malloc(struct dim_light *dim)
{
	dim_light_malloc_xyzl(dim,'x');
	dim_light_malloc_xyzl(dim,'y');
	dim_light_malloc_xyzl(dim,'z');
	dim_light_malloc_xyzl(dim,'l');

}


void dim_light_info_to_buf(struct dat_file *buf,struct dim_light *dim)
{
	long double mul_x=0.0;
	long double mul_y=0.0;
	long double mul_z=0.0;

	get_meter_dim(buf->x_units,&mul_x,dim->x[dim->xlen-1]);
	get_meter_dim(buf->y_units,&mul_y,dim->y[dim->ylen-1]);
	get_meter_dim(buf->z_units,&mul_z,dim->z[dim->zlen-1]);
	buf->y_mul=mul_y;
	buf->x_mul=mul_x;
	buf->z_mul=mul_z;

	strcpy(buf->x_label,_("x-position"));
	strcpy(buf->y_label,_("y-position"));
	strcpy(buf->z_label,_("z-position"));

	buf->x=dim->xlen;
	buf->y=dim->ylen;
	buf->z=dim->zlen;

	buf->logscale_x=0;
	buf->logscale_y=0;

}
