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
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <ray_fun.h>
#include <triangle.h>
#include <triangle_io.h>

/** @file objects.c
	@brief Basic object manipulation
*/

void ray_object_cal_min_max(struct object *obj)
{
	triangles_find_min(&(obj->min),&(obj->tri));
	triangles_find_max(&(obj->max),&(obj->tri));

	//vec_print(&(obj->min));
	//vec_print(&(obj->max));
	//getchar();
}

void ray_object_init(struct object *obj)
{
	//btm
	triangles_init(&(obj->tri));

}

void ray_object_free(struct object *obj)
{
	triangles_free(&(obj->tri));
}

