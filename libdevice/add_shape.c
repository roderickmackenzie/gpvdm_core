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
#include <gpvdm_const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <device.h>
#include <inp.h>
#include <util.h>
#include <triangle_io.h>

/** @file add_shape_to_world.c
	@brief Add a shape to the image
*/

void device_add_shape_to_world(struct simulation *sim,struct device *dev,struct shape *s)
{
	int x=0;
	int z=0;

	double x_pos;
	double y_pos;
	double z_pos;

	struct object *obj;
	struct triangles tri;
	struct vec v;

	char name[200];

	for (x=0;x<s->nx;x++)
	{
		for (z=0;z<s->nz;z++)
		{
			x_pos=s->x0+(s->dx+s->dx_padding)*(double)x;
			y_pos=s->y0;//+s->dy_padding;
			//printf("%Le %s\n",s->y0,s->name);
			//getchar();
			z_pos=s->z0+(s->dz+s->dz_padding)*(double)z;

			triangles_cpy(&tri,&(s->tri));

			triangles_find_min(&v,&tri);

			triangles_sub_vec(&tri,&v);

			triangles_find_max(&v,&tri);

			triangles_div_vec(&tri,&v);

			if (s->flip_y==TRUE)
			{
				v.x=1.0;
				v.y=-1.0;
				v.z=1.0;

				triangles_mul_vec(&tri,&v);

				//v.x=0.0;
				//v.y=1.0;
				//v.z=0.0;

				//triangles_add_vec(&tri,&v);
			}

			v.x=s->dx;
			v.y=s->dy;
			v.z=s->dz;

			triangles_mul_vec(&tri,&v);

			v.x=x_pos;
			v.y=y_pos;
			v.z=z_pos;

			triangles_add_vec(&tri,&v);

			triangles_set_object_type(&tri,RAY_OBJECT);

			obj=ray_add_object(dev,&tri);
			sprintf(name,"%s-(%d,%d)",s->name,x,z);
			strcpy(obj->name,name);

			obj->s=s;

			triangles_free(&tri);

		}
	}

}


