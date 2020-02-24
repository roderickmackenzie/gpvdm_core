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

/** @file dim.h
	@brief Dimention file
*/

#ifndef dim_h
#define dim_h

#include <dat_file_struct.h>

struct dimensions
{
	int zlen;
	int xlen;
	int ylen;


	long double *ymesh;
	long double *xmesh;
	long double *zmesh;

	long double *dy;
	long double *dx;
	long double *dz;

	int srh_bands;

};

struct dim_light
{
	int zlen;
	int xlen;
	int ylen;
	int llen;

	long double *y;
	long double *x;
	long double *z;
	long double *l;

	long double dy;
	long double dx;
	long double dz;
	long double dl;

};

struct dim_heat
{
	int zlen;
	int xlen;
	int ylen;

	long double *y;
	long double *x;
	long double *z;

	long double *dy;
	long double *dx;
	long double *dz;
};

void dim_alloc_gen_untiy_mesh_x(struct dimensions *dim);
void dim_alloc_gen_untiy_mesh_z(struct dimensions *dim);

//dimension
void dim_init(struct dimensions *dim);
void dim_free(struct dimensions *dim);
void dim_alloc(struct dimensions *dim);
void dim_cpy(struct dimensions *out,struct dimensions *in);
void dim_alloc_xyz(struct dimensions *dim,char xyz);
void dim_free_xyz(struct dimensions *dim,char xyz);
void dim_swap(struct dimensions *out,struct dimensions *in);
void dim_info_to_buf(struct dat_file *buf,struct dimensions *dim);

//dim_light
void dim_light_init_xyzl(struct dim_light *dim,char xyzl);
void dim_light_init(struct dim_light *dim);
void dim_light_free_xyzl(struct dim_light *dim,char xyzl);
void dim_light_free(struct dim_light *dim);
void dim_light_malloc_xyzl(struct dim_light *dim,char xyzl);
void dim_light_malloc(struct dim_light *dim);
void dim_light_info_to_buf(struct dat_file *buf,struct dim_light *dim);

//dim_heat
void dim_heat_init_xyz(struct dim_heat *dim,char xyz);
void dim_heat_init(struct dim_heat *dim);
void dim_heat_free_xyz(struct dim_heat *dim,char xyz);
void dim_heat_free(struct dim_heat *dim);
void dim_heat_malloc_xyz(struct dim_heat *dim,char xyz);
void dim_heat_malloc(struct dim_heat *dim);
void dim_heat_info_to_buf(struct dat_file *buf,struct dim_heat *dim);


#endif
