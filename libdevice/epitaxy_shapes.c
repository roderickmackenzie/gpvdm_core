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

/** @file shape.c
	@brief Load the shape files.
*/

#include <epitaxy.h>
#include <sim_struct.h>
#include <shape.h>
#include <string.h>
#include <inp.h>
#include <util.h>
#include <cal_path.h>
#include <i.h>
#include <triangle_io.h>
#include <triangle.h>


void epitaxy_shapes_load(struct simulation *sim,struct epitaxy *in)
{
	int l=0;
	int ii=0;
	char build[100];
	strcpy(build,"");
	int pos=0;
	int len=0;
	struct shape *s;


	for (l=0;l<in->layers;l++)
	{
		in->layer[l].nshape=0;
		if (strcmp(in->shape_file[l],"none")!=0)
		{
			len=strlen(in->shape_file[l])+1;
			for (ii=0;ii<len;ii++)
			{
				if ((in->shape_file[l][ii]==',')||(ii==len-1))
				{
					s=shape_load_file(sim,in,&(in->layer[l].shapes[in->layer[l].nshape]),build);
					in->layer[l].nshape++;
					if (s->flip_y==FALSE)
					{
						s->y0=in->layer[l].y_start;
					}else
					{
						s->y0=in->layer[l].y_stop;
					}

					s->epi_index=l;
					build[0]=0;
					pos=0;
				}else
				{
					build[pos]=in->shape_file[l][ii];
					build[pos+1]=0;
					pos++;
				}

				
			}
		}
		
	}

}


