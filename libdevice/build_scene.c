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
#include <memory.h>
#include <epitaxy_struct.h>
#include <epitaxy.h>
#include <device_fun.h>



/** @file device_build_scene.c
	@brief Build the scene with triangles.
*/


void device_build_scene(struct simulation *sim,struct device *dev)
{

	int l;
	int i;
	int layer;
	int add_layer;
	double xlen=dev->xlen;
	double zlen=dev->zlen;
	int contact_layer=0;
	double dx=xlen*0.01;
	struct epitaxy *epi = &(dev->my_epitaxy);
	//double dz=zlen*0.01;

	double start_z=zlen/2.0;
	double start_x=xlen/2.0;

	double device_height=epitaxy_get_optical_length(epi);
	double sim_window_top=device_height*2.0;
	double scene_y0=device_height*-4.0;
	int c;

	struct shape *s;
	struct object *obj;

	double scene_dx=xlen+dx*2.0;
	double scene_dy=(sim_window_top-scene_y0);
	double scene_dz=dev->zlen;

	s=&(dev->big_box);
	shape_init(sim,s);

	strcpy(s->shape_type,"box");
	strcpy(s->name,"big_box");
	strcpy(s->optical_material,"generic/air");
	s->z0=0.0;
	s->x0=-dx;
	s->y0=scene_y0;

	s->dz=scene_dz;
	s->dx=scene_dx;
	s->dy=scene_dy;
	shape_load_materials(sim,s);
	device_add_shape_to_world(sim,dev,s);
	//obj=add_box(dev,-dx,scene_y0, 0.0,  ,,,RAY_SIM_EDGE);


	//strcpy(obj->name,"big_box");

	//printf("here\n");


	for (l=0;l<epi->layers;l++)
	{

		add_layer=FALSE;
		if (strcmp(epi->layer[l].name,"air")!=0)
		{
			if (epi->layer[l].layer_type==LAYER_CONTACT)
			{
				if (dev->ns.dim.xlen==1)
				{
					add_layer=TRUE;
				}else
				{
					for (c=0;c<dev->ncontacts;c++)
					{
						if ((dev->contacts[c].position==TOP)&&(l==0))
						{
							//printf("a> %d %s %d\n",l,dev->contacts[c].name,dev->contacts[c].position);
							//getchar();
							dev->contacts[c].shape.dy=epi->layer[l].width;
							dev->contacts[c].shape.dz=dev->zlen;
							dev->contacts[c].shape.y0=epi->layer[l].y_start;
							//dev->contacts[c].shape.x0=0.0;
							dev->contacts[c].shape.z0=0.0;
							dev->contacts[c].shape.nx=1;
							dev->contacts[c].shape.nz=1;
							dev->contacts[c].shape.flip_y=FALSE;
							device_add_shape_to_world(sim,dev,&(dev->contacts[c].shape));		//,epi->layer[l].y_stop
						}

						if ((dev->contacts[c].position==BOTTOM)&&(l==epi->layers-1))
						{
//							printf("b> %d %s %d\n",l,dev->contacts[c].name,dev->contacts[c].position);
							dev->contacts[c].shape.dy=epi->layer[l].width;
							dev->contacts[c].shape.dz=dev->zlen;
							dev->contacts[c].shape.z0=0.0;
							dev->contacts[c].shape.y0=epi->layer[l].y_start;
							dev->contacts[c].shape.x0=0.0;
							dev->contacts[c].shape.nx=1;
							dev->contacts[c].shape.nz=1;
							dev->contacts[c].shape.flip_y=FALSE;
							device_add_shape_to_world(sim,dev,&(dev->contacts[c].shape));	//,epi->layer[l].y_stop

						}


						//

					}
				}
				contact_layer++;
			}else
			{
				add_layer=TRUE;
			}

			if (add_layer==TRUE)
			{
				//printf("%Le\n",epi->layer[l].y_start);
				//getchar();
				obj=add_box(dev,0.0,epi->layer[l].y_start,0.0,xlen,fabs(epi->layer[l].width),zlen,RAY_OBJECT);
				obj->epi_layer=l;
				strcpy(obj->name,epi->layer[l].name);
			}
		}
	}


	for (c=0;c<dev->ncontacts;c++)
	{
		if ((dev->contacts[c].position==LEFT))
		{
		//	printf("a> %d %s %d\n",l,dev->contacts[c].name,dev->contacts[c].position);
			//dev->contacts[c].shape.dy=epi->layer[l].width;
			//dev->contacts[c].shape.z0=0.0;
			//dev->contacts[c].shape.x0=1e-5;
			//dev->contacts[c].shape.nx=1;
			//dev->contacts[c].shape.nz=1;

			device_add_shape_to_world(sim,dev,&(dev->contacts[c].shape));
		}
	}


	for (l=0;l<epi->layers;l++)
	{
		for (i=0;i<epi->layer[l].nshape;i++)
		{
			s=&epi->layer[l].shapes[i];
			s->x0=0.0;
			device_add_shape_to_world(sim,dev,s);
		}
	}


	device_dump_world_to_file(sim,dev,"device.dat");
}
