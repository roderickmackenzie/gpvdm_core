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
#include "const.h"
#include "light.h"
#include "device.h"
#include "const.h"
#include "dump.h"
#include "config.h"
#include "inp.h"
#include "util.h"
#include "hard_limit.h"
#include "lang.h"
#include "log.h"
#include "memory.h"


static int unused __attribute__((unused));


void light_init_mesh(struct simulation *sim,struct light *in,struct epitaxy *my_epitaxy)
{
printf_log(sim,"init: mesh\n");
	int i;
	gdouble pos=0.0;
	pos=in->dx;

	int layer=0;
	gdouble layer_end=in->thick[layer];

	for (i=0;i<in->points;i++)
	{
		in->x[i]=pos;
		in->layer_end[i]=layer_end-pos;
		in->layer[i]=layer;
		if (in->device_start_layer>=layer) in->device_start_i=i;

		pos+=in->dx;

		if (pos>layer_end)
		{

			//do
			//{
			if (layer<(in->layers-1))
			{
				layer++;
			//}while(in->thick[layer]==0.0);

				layer_end=layer_end+in->thick[layer];
			}
		}

	}
	in->device_start_i++;

	in->dl=(in->lstop-in->lstart)/((gdouble)in->lpoints);

	pos=in->lstart;
	for (i=0;i<in->lpoints;i++)
	{
		in->l[i]=pos;
		pos+=in->dl;
	}

	int ii;

	for (i=0;i<in->lpoints;i++)
	{
		for (ii=0;ii<in->points;ii++)
		{

			in->alpha[i][ii]=inter_get_noend(&(my_epitaxy->mat[in->layer[ii]]),in->l[i]);
			in->alpha0[i][ii]=in->alpha[i][ii];
			in->n[i][ii]=inter_get_noend(&(my_epitaxy->mat_n[in->layer[ii]]),in->l[i]);

		}
	}

	if (in->flip_field==TRUE)
	{
		for (i=0;i<in->lpoints;i++)
		{
			memory_flip_1d_long_double(in->alpha[i],in->points);
			memory_flip_1d_long_double(in->alpha0[i],in->points);
			memory_flip_1d_long_double(in->n[i],in->points);
		}
	}



	light_calculate_complex_n(in);


	for (i=0;i<in->lpoints;i++)
	{
		in->sun_norm[i]=inter_get_hard(&(in->sun_read),in->l[i]);
	}

	gdouble tot=0.0;
	for  (i=0;i<in->lpoints;i++)
	{
		tot+=in->dl*in->sun_norm[i];
	}

	//for  (i=0;i<in->lpoints;i++)
	//{
	//	in->sun_norm[i]/=tot;
	//}


}
