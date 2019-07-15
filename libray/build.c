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
#include <const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>
#include <device.h>
#include <inp.h>
#include <util.h>

/** @file build.c
	@brief Set up the simulation window for the ray tracer
*/

//in->l[lam]
void light_update_ray_mat(struct simulation *sim,struct epitaxy *my_epitaxy,struct image *my_image,long double lambda)
{
	int i;
	int layer;
	long double n=1.0;
	long double alpha=0.0;

	for (i=0;i<my_image->objects;i++)
	{
		layer=my_image->obj_mat_number[i];

		if (layer==-1)
		{
			alpha=1e-3;
			n=1.0;
		}else
		{
			alpha=inter_get_noend(&(my_epitaxy->mat[layer]),lambda);
			n=inter_get_noend(&(my_epitaxy->mat_n[layer]),lambda);
		}

		my_image->obj_n[i]=n;
		my_image->obj_alpha[i]=alpha;
		//printf("%Le %Le %d %d\n",n,alpha,i,layer);
		//getchar();
	}



}


void ray_read_config(struct simulation *sim,struct image *my_image)
{
	int i;
	struct inp_file inp;
	char temp[200];

	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_input_path(sim),"ray.inp");

	inp_check(sim,&inp,1.0);

	inp_search_int(sim,&inp,&(my_image->theta_steps),"#ray_theta_steps");
	inp_search_int(sim,&inp,&(my_image->ray_wavelength_points),"#ray_wavelength_points");

	inp_search_int(sim,&inp,&(my_image->escape_bins),"#ray_escape_bins");

	inp_search_gdouble(sim,&inp,&(my_image->ray_xsrc),"#ray_xsrc");
	inp_search_gdouble(sim,&inp,&(my_image->ray_ysrc),"#ray_ysrc");
	inp_search_gdouble(sim,&inp,&(my_image->ray_zsrc),"#ray_zsrc");

	inp_search_gdouble(sim,&inp,&(my_image->ray_theta_start),"#ray_theta_start");
	inp_search_gdouble(sim,&inp,&(my_image->ray_theta_stop),"#ray_theta_stop");

	inp_search_string(sim,&inp,temp,"#ray_auto_run");
	my_image->ray_auto_run=english_to_bin(sim,temp);

	inp_search_string(sim,&inp,temp,"#ray_input_spectrum");

	join_path(3, my_image->input_spectrum_file, sim->emission_path, temp,"spectra.inp");

	if (inp_isfile(sim,my_image->input_spectrum_file)!=0)
	{
		ewe(sim,"The emission file %s does not exist",my_image->input_spectrum_file);
	}

	ray_load_emission(sim,my_image);

	inp_free(sim,&inp);

	my_image->lam=malloc(sizeof(long double)*my_image->ray_wavelength_points);
	my_image->extract_eff=malloc(sizeof(long double)*my_image->ray_wavelength_points);

	my_image->ang_escape=(long double **)malloc(sizeof(long double*)*my_image->ray_wavelength_points);
	my_image->angle=(long double *)malloc(sizeof(long double)*my_image->escape_bins);

	long double da=180.0/(long double)my_image->escape_bins;
	long double apos=0.0;
	for (i=0;i<my_image->escape_bins;i++)
	{
		apos+=da;
		my_image->angle[i]=apos;
	}

	long double lam=my_image->input_spectrum.x[0];
	long double dl=(my_image->input_spectrum.x[my_image->input_spectrum.len-1]-my_image->input_spectrum.x[0])/((long double)my_image->ray_wavelength_points);

	for (i=0;i<my_image->ray_wavelength_points;i++)
	{
		my_image->lam[i]=lam;
		my_image->ang_escape[i]=(long double*)malloc(sizeof(long double)*my_image->escape_bins);

		lam+=dl;
	}

}

void ray_free(struct simulation *sim,struct image *my_image)
{
	int i=0;

	if (my_image->lam!=NULL)
	{
		free(my_image->lam);
	}

	if (my_image->extract_eff!=NULL)
	{
		free(my_image->extract_eff);
	}

	for (i=0;i<my_image->ray_wavelength_points;i++)
	{
		free(my_image->ang_escape[i]);
	}

	free(my_image->ang_escape);
	free(my_image->angle);

	inter_free(&(my_image->input_spectrum));
}

void light_setup_ray(struct simulation *sim,struct device *cell,struct image *my_image,struct epitaxy *my_epitaxy)
{

	int i;

	double xlen=cell->xlen;
	double dx=xlen*0.01;
	double device_start=epitaxy_get_device_start(my_epitaxy);
	double device_stop=epitaxy_get_device_stop(my_epitaxy);

	double start_y=(device_stop-device_start)/2.0;

	double device_height=epitaxy_get_optical_length(my_epitaxy);
	double sim_window_top=device_height*4.0;
	double sim_window_btm=device_height*-2.0;

	my_image->y_escape_level=-device_height*1.1;

	add_box(my_image,0.0,sim_window_btm,xlen+dx*2.0,sim_window_top,-1,TRUE);

	for (i=0;i<my_epitaxy->layers;i++)
	{
		add_box(my_image,dx,my_epitaxy->y_pos[i],xlen,fabs(my_epitaxy->width[i]),i,FALSE);
	}


	my_image->n_start_rays=10;
	double x_start=dx+dx/2.0;
	double x_stop=dx+xlen-dx/2.0;
	dx=(x_stop-x_start)/((double)my_image->n_start_rays);
	double x_pos=x_start;
	
	if (my_image->ray_xsrc==-1.0)
	{
		for (i=0;i<my_image->n_start_rays;i++)
		{
			my_image->start_rays[i].x=x_pos;
			my_image->start_rays[i].y=start_y;
			my_image->start_rays[i].z=0.0;

			x_pos=x_pos+dx;
		}
	}else
	{
		for (i=0;i<my_image->n_start_rays;i++)
		{
			my_image->start_rays[i].x=my_image->ray_xsrc;
			my_image->start_rays[i].y=my_image->ray_ysrc;
			my_image->start_rays[i].z=my_image->ray_zsrc;
			x_pos=x_pos+dx;
		}
	}

	//dump_plane(my_image);
	//dump_plane_to_file(my_image);
}
