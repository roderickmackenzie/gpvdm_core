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

/** @file ray_search_intersect.c
	@brief Ray tracing for the optical model, this should really be split out into it's own library.
*/
int search_triangle(struct simulation *sim,struct image *in,struct ray *my_ray)
{
int i;
int found=-1;
int pos=-1;

double min_dist=1000.0;
double dist=0;

struct vec ret;
vec_init(&ret);

struct vec tmp;
vec_init(&tmp);

struct vec store;
vec_init(&store);

	for (i=0;i<in->lines;i++)
	{
		//dump_plane(sim,in);
		//dump_plane_to_file("lines.dat",in);
		found=ray_intersect(&ret,&(in->p[i]),my_ray);
		//vec_print(&ret);
		//printf("found=%d\n",found);
		//getchar();

		if (found==TRUE)
		{
			vec_cpy(&tmp,&ret);
			//vec_print(&ret);
			//getchar();
			//vec_print(&(in->rays[ray].xy));			
			vec_sub(&tmp,&(my_ray->xy));
			dist=vec_fabs(&tmp);

			if (dist>1e-12)
			{
				if (dist<min_dist)
				{
					pos=i;
					vec_cpy(&(my_ray->xy_end),&ret);
					min_dist=dist;
				}
			}
		}
		
	}

//dump_plane_to_file("lines.dat",in);
//printf("exit\n");
//getchar();
return pos;	
}

int search_object(struct simulation *sim,struct image *in,struct ray *my_ray)
{
	int t=0;
	t=search_triangle(sim,in,my_ray);
	if (t==-1)
	{
		return t;
	}

	return in->p[t].id;
}
