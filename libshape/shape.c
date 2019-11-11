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


int shape_get_index(struct simulation *sim,struct epitaxy *in,long double x,long double y,long double z)
{
	int i=0;
	int l=0;
	struct shape *s;

	for (l=0;l<in->layers;l++)
	{
		for (i=0;i<in->layer[l].nshape;i++)
		{
			s=&in->layer[l].shapes[i];

			if ((x>s->x0)&&(x<(s->x0+s->dx)))
			{
				if ((y>s->y0)&&(y<((s->y0)+(s->dy))))
				{
					//printf("%Le %Le %Le \n",z,s->z0,s->dz);

					//if ((z>s->z0)&&(z<((s->z0)+(s->dz))))
					{
						return s->dos_index;
					}
				}

			}
		}
	}
	return -1;
	
}

void shape_free(struct simulation *sim,struct shape *s)
{
	if (strcmp(s->optical_material,"none")!=0)
	{
		inter_free(&(s->n));
		inter_free(&(s->alpha));
	}

	triangles_free((&(s->tri)));
}

void shape_free_materials(struct simulation *sim,struct epitaxy *in)
{
	int i;
	int l;
	struct shape *s;

	for (l=0;l<in->layers;l++)
	{
		for (i=0;i<in->layer[l].nshape;i++)
		{
			s=&in->layer[l].shapes[i];
			shape_free(sim,s);
		}
	}
}

struct shape *shape_load_file(struct simulation *sim,struct epitaxy *in,struct shape *s, char *file_name)
{
	char full_file_name[PATH_MAX];
	char file_path[PATH_MAX];
	char temp_path[PATH_MAX];

	sprintf(full_file_name,"%s.inp",file_name);

	printf("Loading shape file %s\n",full_file_name);

	struct inp_file inp;

	inp_init(sim,&inp);

	s->dos_index=0;

	if (inp_load(sim,&inp,full_file_name)==0)
	{
		inp_search_gdouble(sim,&inp,&(s->dx),"#shape_dx");
		inp_search_gdouble(sim,&inp,&(s->dy),"#shape_dy");
		inp_search_gdouble(sim,&inp,&(s->dz),"#shape_dz");

		inp_search_gdouble(sim,&inp,&(s->dx_padding),"#shape_padding_dx");
		inp_search_gdouble(sim,&inp,&(s->dy_padding),"#shape_padding_dy");
		inp_search_gdouble(sim,&inp,&(s->dz_padding),"#shape_padding_dz");


		inp_search_int(sim,&inp,&(s->nx),"#shape_nx");
		inp_search_int(sim,&inp,&(s->nz),"#shape_nz");

		inp_search_string(sim,&inp,s->name,"#shape_name");
		inp_search_string(sim,&inp,s->shape_type,"#shape_type");

		join_path(3,temp_path,get_shape_path(sim),s->shape_type,"shape.inp");
		triangle_load_from_file(sim,(&s->tri),temp_path);

		inp_search_string(sim,&inp,s->optical_material,"#shape_optical_material");

		inp_search_string(sim,&inp,s->dos_file,"#shape_dos");

		inp_search_gdouble(sim,&inp,&(s->x0),"#shape_x0");
		inp_search_gdouble(sim,&inp,&(s->z0),"#shape_z0");

		if (strcmp(s->dos_file,"none")!=0)
		{
			s->dos_index=in->electrical_layers;
			epitaxy_load_dos_files(sim,in, s->dos_file,"none","none");
		}

		if (strcmp(s->optical_material,"none")!=0)
		{
			join_path(3, file_path,get_materials_path(sim),s->optical_material,"alpha.gmat");
			if (isfile(file_path)!=0)
			{
				ewe(sim,"File %s not found",file_path);
			}

			inter_load(sim,&(s->alpha),file_path);
			inter_sort(&(s->alpha));

			join_path(3, file_path,get_materials_path(sim),s->optical_material,"n.gmat");
			if (isfile(file_path)!=0)
			{
				ewe(sim,"File %s not found",file_path);
			}

			inter_load(sim,&(s->n),file_path);
			inter_sort(&(s->n));

		}

		inp_free(sim,&inp);
	}else
	{
		ewe(sim,"file %s not found\n",full_file_name);
	}

return s;
}

void shape_load(struct simulation *sim,struct epitaxy *in)
{
	int i=0;
	int ii=0;
	char build[100];
	strcpy(build,"");
	int pos=0;
	int len=0;
	struct shape *s;


	for (i=0;i<in->layers;i++)
	{
		in->layer[i].nshape=0;
		if (strcmp(in->shape_file[i],"none")!=0)
		{
			len=strlen(in->shape_file[i])+1;
			for (ii=0;ii<len;ii++)
			{
				if ((in->shape_file[i][ii]==',')||(ii==len-1))
				{
					s=shape_load_file(sim,in,&(in->layer[i].shapes[in->layer[i].nshape]),build);
					in->layer[i].nshape++;
					s->y0=in->y_pos[i];
					s->epi_index=i;
					build[0]=0;
					pos=0;
				}else
				{
					build[pos]=in->shape_file[i][ii];
					build[pos+1]=0;
					pos++;
				}

				
			}
		}
		
	}
}


