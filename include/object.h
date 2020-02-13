// 
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
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
// 

/** @file object.h
@brief ray tracing header files.
*/
#ifndef object_h
#define object_h

#include <vec.h>
#include <sim_struct.h>
#include <triangle.h>

struct object
{
	int epi_layer;
	char name[100];
	int uid;
	double *n;
	double *alpha;
	struct triangles tri;
	struct vec min;
	struct vec max;
	struct shape* s;		//This is a poinnter to the origonal shape which generated the object
};

//Object
void object_flip_y_axis(struct object *obj);
void object_sub_y(struct object *obj,double y);
void object_add_y(struct object *obj,double y);
double object_get_min_y(struct object *obj);
void object_init(struct object *obj);
void object_free(struct object *obj);
void object_cal_min_max(struct object *obj);
void object_nalpha_malloc(struct object *obj,int l_max);
void object_nalpha_free(struct object *obj);
void object_malloc(struct object *obj);
#endif
