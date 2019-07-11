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
#include <sim_struct.h>

#define WAIT 0
#define READY 1
#define DONE 2

#define TRUE 1
#define FALSE 0

#define RAY_MAX 5000

struct plane
{
	struct vec xy0;
	struct vec xy1;
	int id;
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

struct image
{
	struct plane p[1000];
	struct ray rays[RAY_MAX];
	int lines;
	int nrays;
	int obj_mat_number[100];
	double obj_n[100];
	double obj_alpha[100];
	int objects;
	struct vec start_rays[100];
	int n_start_rays;
	int theta_steps;
	double y_escape_level;

	int escape_angle_bins;
	long double *angle;
	long double **ang_escape;
	int ray_wavelength_points;
	long double ray_wavelength_start;
	long double ray_wavelength_stop;
	long double *lam;
	long double *extract_eff;
	int ray_auto_run;
	int ray_escape_bins;
};


#endif
