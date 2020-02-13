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

/** @file light_mesh.c
	@brief Performs meshing for the light model.
*/

#include "util.h"
#include "gpvdm_const.h"
#include "light.h"
#include "device.h"
#include "gpvdm_const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "lang.h"
#include "log.h"
#include "memory.h"
#include <light_fun.h>
#include <epitaxy.h>
#include <ray_fun.h>

static int unused __attribute__((unused));


void light_build_mesh(struct simulation *sim,struct light *li,struct dimensions *electrical_dim)
{
	int y=0;
	int x=0;
	int z=0;
	int l=0;

	long double pos=0.0;

	struct dim_light *dim=&(li->dim);

	printf_log(sim,"init: mesh\n");
	for (x=0;x<electrical_dim->xlen;x++)
	{
		dim->x[x]=electrical_dim->xmesh[x];
	}

	for (z=0;z<electrical_dim->zlen;z++)
	{
		dim->z[z]=electrical_dim->zmesh[z];
	}

	pos=0.0;
	for (y=0;y<dim->ylen;y++)
	{
		dim->y[y]=pos;
		pos=pos+dim->dy;
		//printf("%Le %Le\n",pos,dim->dy);
		//getchar();

	}

	pos=li->lstart;
	for (l=0;l<dim->llen;l++)
	{
		dim->l[l]=pos;
		pos+=dim->dl;
	}

}


void light_build_obj_pointer_array(struct simulation *sim,struct light *li, struct device *dev)
{
	printf_log(sim,"Light: Building object pointer array\n");
	int x=0;
	int y=0;
	int z=0;

	struct dim_light *dim=&(li->dim);
	struct vec v;
	//for (x=0;x<dim->xlen;x++)
	//{
	//	printf("%d %Le\n",x,dim->x[x]);
	//}
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				v.z=dim->z[z];
				v.x=dim->x[x];
				v.y=dim->y[y];

				li->obj[z][x][y]=ray_obj_search_xyz(sim,dev,&v);
				//printf("%Le %Le %Le %d\n",dim->z[z],dim->x[x],dim->y[y],li->obj[z][x][y]);
				//if (li->obj[z][x][y]!=0)
				//{
				//	printf("%s\n",li->obj[z][x][y]->name);
				//}
				//getchar();

			}
		}
	}


}

void light_build_materials_arrays(struct simulation *sim,struct light *li, struct device *dev)
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;
	int layer=0;
	//long double ypos=0.0;
	long double n=0.0;
	long double alpha=0.0;
	long double lam=0.0;
	struct object *obj;
	struct epitaxy *epi=li->epi;
	struct dim_light *dim=&(li->dim);
	struct shape *s;
	//struct vec v;

	light_build_obj_pointer_array(sim,li, dev);

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				//v.z=dim->z[z];
				//v.x=dim->x[x];
				//v.y=dim->y[y];

				obj=li->obj[z][x][y];//ray_obj_search_xyz(sim,dev,&v);
				//if (obj->epi_layer==-1)
				//{
				//	printf("%s\n",obj->s->optical_material);
				//}

				//printf("%s\n",li->obj[z][x][y]->name);
				//getchar();

				for (l=0;l<dim->llen;l++)
				{
					lam=dim->l[l];


					if (obj->epi_layer!=-1)
					{
						layer=obj->epi_layer;//epitaxy_get_layer(epi,ypos);
						n=inter_get_noend(&(epi->layer[layer].n),lam);
						alpha=inter_get_noend(&(epi->layer[layer].alpha),lam);

						li->alpha[z][x][y][l]=alpha;
						li->alpha0[z][x][y][l]=alpha;
						li->n[z][x][y][l]=n;
					}else
					{
						s=obj->s;

						if (strcmp(s->optical_material,"none")==0)
						{
							ewe(sim,"No optical material defined for object: %s",obj->name);
						}
						//if (li->obj[z][x][y]!=0)
						//{

						//}

						n=inter_get_noend(&(s->n),lam);
						alpha=inter_get_noend(&(s->alpha),lam);
						li->alpha[z][x][y][l]=alpha;
						li->alpha0[z][x][y][l]=alpha;
						li->n[z][x][y][l]=n;

					}
				}
			}
			//getchar();
		}
	}

		//getchar();
	if (li->flip_field==TRUE)
	{
		flip_light_zxyl_long_double_y(sim, dim,li->alpha);
		flip_light_zxyl_long_double_y(sim, dim,li->alpha0);
		flip_light_zxyl_long_double_y(sim, dim,li->n);
	}


	light_calculate_complex_n(li);


	for (l=0;l<dim->llen;l++)
	{
		//printf("%Le\n",dim->l[l]);
		//getchar();
		li->sun_norm[l]=inter_get_hard(&(li->sun_read),dim->l[l]);
	}


}
