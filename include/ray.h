// 
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
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
// 

/** @file ray.h
@brief ray tracing header files.
*/
#ifndef ray_h
#define ray_h

#include <vec.h>
#include <i.h>
#include <sim_struct.h>

#define WAIT 0
#define READY 1
#define DONE 2

#define TRUE 1
#define FALSE 0

struct triangle
{
	struct vec xy0;
	struct vec xy1;
	struct vec xy2;
	int object_uid;
	int edge;
};

struct ray
{
	struct vec xy;
	struct vec xy_end;
	struct vec dir;
	int state;
	int bounce;
	double mag;
};

struct object
{
	int epi_layer;
	int shape_number;
	char name[100];
	double n;
	int uid;
	double alpha;
};

struct image
{
	int start_of_shapes;
	int triangles_max;
	int triangles;
	struct triangle *tri;
	struct ray *rays;
	int nrays;
	int nray_max;
	struct object obj[1000];

	int objects;
	struct vec start_rays[100];
	int n_start_rays;
	int top_of_done_rays;

	double y_escape_level;
	char input_spectrum_file[PATH_MAX];
	struct istruct input_spectrum;
	long double *angle;
	long double **ang_escape;
	int ray_wavelength_points;
	double *lam;
	double cur_lam;
	double *extract_eff;
	double avg_extract_eff;
	int ray_auto_run;
	int escape_bins;
	double ray_xsrc;
	double ray_ysrc;
	double ray_zsrc;

	int theta_steps;
	double ray_theta_start;
	double ray_theta_stop;


	int phi_steps;
	double ray_phi_start;
	double ray_phi_stop;

};


#endif
