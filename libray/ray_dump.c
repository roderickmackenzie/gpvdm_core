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
#include <buffer.h>

/** @file ray.c
	@brief Ray tracing for the optical model, this should really be split out into it's own library.
*/

void dump_plane_to_file(char *file_name,struct image *in)
{
	FILE *out;

	out=fopen("lines.dat","w");
	int i=0;

	for (i=0;i<in->lines;i++)
	{
		fprintf(out,"%le %le\n",in->p[i].xy0.x,in->p[i].xy0.y);
		fprintf(out,"%le %le\n",in->p[i].xy1.x,in->p[i].xy1.y);
		fprintf(out,"\n");
		fprintf(out,"\n");

	}

	fclose(out);
	//file_name
	out=fopen(file_name,"a");

	for (i=0;i<in->nrays;i++)
	{
		if (in->rays[i].state==DONE)
		{
			fprintf(out,"%le %le\n",in->rays[i].xy.x,in->rays[i].xy.y);
			fprintf(out,"%le %le\n",in->rays[i].xy_end.x,in->rays[i].xy_end.y);
			fprintf(out,"\n");
		}
		
	}

	fclose(out);
	
	//out=fopen("start.out","w");
	//for (i=0;i<in->n_start_rays;i++)
	//{
	//	fprintf(out,"%le %le\n\n",in->start_rays[i].x,in->start_rays[i].y);
	//}
	//fclose(out);
	
}

void dump_plane(struct simulation *sim,struct image *in)
{
	int i=0;
	printf_log(sim,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	for (i=0;i<in->n_start_rays;i++)
	{
		printf_log(sim,"%le %le\n",in->start_rays[i].x,in->start_rays[i].y);
	}

	printf_log(sim,"lines:\n");
	for (i=0;i<in->lines;i++)
	{
		printf_log(sim,"%le %le %le %le %d\n",in->p[i].xy0.x,in->p[i].xy0.y,in->p[i].xy1.x,in->p[i].xy1.y,in->p[i].edge);


	}

	printf_log(sim,"rays x,y,x_vec,y_vec:\n");


	for (i=0;i<in->nrays;i++)
	{
		printf_log(sim,"%d (%le,%le) (%le,%le) %lf %lf mag=%lf\n",in->rays[i].state,in->rays[i].xy.x,in->rays[i].xy.y,in->rays[i].xy_end.x,in->rays[i].xy_end.y,in->rays[i].dir.x,in->rays[i].dir.y,in->rays[i].mag);
		
	}
	printf_log(sim,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

}

void dump_extraction_efficiency(struct simulation *sim,struct image *in)
{
	int i;
	char temp[200];
	struct buffer buf;
	buffer_init(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1.0;
	buf.y_mul=1e9;
	strcpy(buf.title,"Photon escape probability");
	strcpy(buf.type,"linegraph");
	strcpy(buf.y_label,"Wavelength");
	strcpy(buf.data_label,"Probability");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"a.u.");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=1;
	buf.y=in->ray_wavelength_points;
	buf.z=1;
	buffer_add_info(sim,&buf);

	for (i=0;i<in->ray_wavelength_points;i++)
	{
		sprintf(temp,"%Le %Le\n",in->lam[i],in->extract_eff[i]);
		buffer_add_string(&buf,temp);
	}

	buffer_dump_path(sim,"","escape_probability.dat",&buf);
	buffer_free(&buf);
}

void dump_ang_escape(struct simulation *sim,struct image *in)
{
	int x;
	int y;
	char temp[200];
	struct buffer buf;
	buffer_init(&buf);

	buffer_malloc(&buf);
	buf.y_mul=1.0;
	buf.y_mul=1e9;
	strcpy(buf.title,"Photon escape probability");
	strcpy(buf.type,"heat");
	strcpy(buf.y_label,"Angle");
	strcpy(buf.x_label,"Wavelength");

	strcpy(buf.data_label,"Probability");
	strcpy(buf.y_units,"Degrees");
	strcpy(buf.y_units,"nm");
	strcpy(buf.data_units,"a.u.");
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->escape_angle_bins;
	buf.y=in->ray_wavelength_points;
	buf.z=1;
	buffer_add_info(sim,&buf);

	for (y=0;y<in->escape_angle_bins;y++)
	{
		for (x=0;x<in->ray_wavelength_points;x++)
		{
			sprintf(temp,"%Le %Le %Le\n",in->lam[x],my_image->angle[y],in->ang_escape[x][y]);
			buffer_add_string(&buf,temp);
		}
		buffer_add_string(&buf,"\n");
	}

	buffer_dump_path(sim,"","ang_escape.dat",&buf);
	buffer_free(&buf);
}
