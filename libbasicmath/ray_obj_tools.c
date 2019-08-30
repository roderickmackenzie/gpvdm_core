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

#include <stdio.h>
#include <ray.h>
#include <ray_fun.h>
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <device.h>
#include <inp.h>
#include <util.h>

/** @file ray_obj_tools.c
	@brief Object tools
*/

double ray_tri_get_min_y(struct triangle* tri)
{
	double min=tri->xy0.y;
	if (min>tri->xy1.y)
	{
		min=tri->xy1.y;
	}

	if (min>tri->xy2.y)
	{
		min=tri->xy2.y;
	}

	return min;
}

void ray_tri_flip_y_axis(struct triangle* tri,double y)
{
	vec_flip_y_axis(&(tri->xy0),y);
	vec_flip_y_axis(&(tri->xy1),y);
	vec_flip_y_axis(&(tri->xy2),y);
}

void ray_tri_sub_y(struct triangle* tri,double y)
{
	struct vec a;
	vec_set(&a,0,y,0);
	vec_sub(&(tri->xy0),&a);
	vec_sub(&(tri->xy1),&a);
	vec_sub(&(tri->xy2),&a);

}

void ray_tri_add_y(struct triangle* tri,double y)
{
	struct vec a;
	vec_set(&a,0,y,0);
	vec_add(&(tri->xy0),&a);
	vec_add(&(tri->xy1),&a);
	vec_add(&(tri->xy2),&a);

}

double ray_obj_get_min_y(struct simulation *sim,struct image *in,struct object *obj)
{
	int i=0;
	double min=1e9;//ray_tri_get_min_y(&(in->tri[0]));
	double min_new=min;
	struct triangle *tri;
	for (i=0;i<obj->tri.len;i++)
	{
		tri=&(obj->tri.data[i]);
		min_new=ray_tri_get_min_y(tri);
		if (min_new<min)
		{
			min=min_new;
		}
	}

return min;
}


void ray_object_flip_y_axis(struct simulation *sim,struct image *in,struct object *obj)
{
	int i=0;
	double y=ray_obj_get_min_y(sim,in,obj);
	struct triangle *tri;

	for (i=0;i<obj->tri.len;i++)
	{
		tri=&(obj->tri.data[i]);
		ray_tri_flip_y_axis(tri,y);
	}

}

void ray_object_sub_y(struct simulation *sim,struct image *in,struct object *obj,double y)
{
	int i=0;
	struct triangle *tri;
	for (i=0;i<obj->tri.len;i++)
	{
		tri=&(obj->tri.data[i]);
		ray_tri_sub_y(tri,y);
	}

}

void ray_object_add_y(struct simulation *sim,struct image *in,struct object *obj,double y)
{
	int i=0;
	struct triangle *tri;
	for (i=0;i<obj->tri.len;i++)
	{
		tri=&(obj->tri.data[i]);
		ray_tri_add_y(tri,y);
	}

}
