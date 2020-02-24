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
#include <triangles.h>
#include <memory.h>
#include <epitaxy_struct.h>
#include <epitaxy.h>
#include <dat_file.h>

/** @file scene_dump.c
	@brief Dump the scene to file
*/


void device_dump_world_to_file(struct simulation *sim,struct device *dev,char *file_name)
{
	int i;

	char temp[200];

	//printf("file dump\n");
	struct dat_file buf;
	buffer_init(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	strcpy(buf.title,"Ray trace triange file");
	strcpy(buf.type,"poly");
	strcpy(buf.y_label,"Position");
	strcpy(buf.x_label,"Position");
	strcpy(buf.data_label,"Position");

	strcpy(buf.y_units,"m");
	strcpy(buf.x_units,"m");
	strcpy(buf.data_units,"m");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=1;
	buf.y=dev->triangles;
	buf.z=1;
	buffer_add_info(sim,&buf);
	struct object *obj;
	struct triangle *tri;


	int o=0;
	for (o=0;o<dev->objects;o++)
	{
		obj=&(dev->obj[o]);

		sprintf(temp,"#name %s\n",obj->name);
		buffer_add_string(&buf,temp);

		for (i=0;i<obj->tri.len;i++)
		{
			tri=&(obj->tri.data[i]);
			sprintf(temp,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
			buffer_add_string(&buf,temp);

			sprintf(temp,"%le %le %le\n",tri->xy1.z,tri->xy1.x,tri->xy1.y);
			buffer_add_string(&buf,temp);

			sprintf(temp,"%le %le %le\n",tri->xy2.z,tri->xy2.x,tri->xy2.y);
			buffer_add_string(&buf,temp);

			sprintf(temp,"%le %le %le\n",tri->xy0.z,tri->xy0.x,tri->xy0.y);
			buffer_add_string(&buf,temp);

			sprintf(temp,"\n");
			buffer_add_string(&buf,temp);

			sprintf(temp,"\n");
			buffer_add_string(&buf,temp);


		}

	}

	buffer_dump_path(sim,"",file_name,&buf);
	buffer_free(&buf);
}

