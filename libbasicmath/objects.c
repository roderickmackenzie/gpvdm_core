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
#include <gpvdm_const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <ray_fun.h>
#include <util.h>
#include <triangle.h>
#include <triangles.h>


/** @file ray_search_intersect.c
	@brief Ray tracing for the optical model, this should really be split out into it's own library.
*/

void objects_dump(struct simulation *sim,struct device *dev)
{
int o=0;

struct object *obj;

	for (o=0;o<dev->objects;o++)
	{
		obj=&(dev->obj[o]);
		printf("%d:%s %d",o,obj->name,obj->tri.len);
		printf("\t(%le,%le,%le)",obj->min.x,obj->min.y,obj->min.z);
		printf("\t(%le,%le,%le)\n",obj->max.x,obj->max.y,obj->max.z);
	}

}

