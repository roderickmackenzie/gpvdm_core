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


/** @file epitaxy.h
	@brief Read the epitaxy from the epitaxy.inp file.
*/


#ifndef epitaxy_h
#define epitaxy_h
#include "advmath.h"
#include <sim_struct.h>
#include <shape_struct.h>
#include <epitaxy_struct.h>
#include <device.h>

void epitaxy_load(struct simulation *sim,struct epitaxy *in, char *file);
gdouble epitaxy_get_electrical_length(struct epitaxy *in);
gdouble epitaxy_get_optical_length(struct epitaxy *in);
int epitaxy_get_layer(struct epitaxy *in,gdouble pos);
int epitaxy_get_electrical_material_layer(struct epitaxy *in,gdouble pos);
gdouble epitaxy_get_device_start(struct epitaxy *in);
gdouble epitaxy_get_device_stop(struct epitaxy *in);
gdouble epitaxy_get_device_start_i(struct epitaxy *in);
int epitaxy_get_epitaxy_layer_using_electrical_pos(struct epitaxy *in,gdouble pos);
void epitaxy_load_materials(struct simulation *sim,struct epitaxy *in);
void epitaxy_free(struct simulation *sim,struct epitaxy *in);
void epitaxy_free_materials(struct epitaxy *in);
void epitaxy_load_dos_files(struct simulation *sim,struct epitaxy *in, char *dos_file,char *lumo_file,char *homo_file);
void epitaxy_load_emission(struct simulation *sim,struct epi_layer *layer);
void epitaxy_load_electrical_file(struct simulation *sim,struct epi_layer *layer);
void epitaxy_mask(struct simulation *sim,struct device *dev);
void epitaxy_shapes_load(struct simulation *sim,struct epitaxy *in);

//optical
long double epitaxy_get_optical_problem_start(struct epitaxy *in);
long double epitaxy_get_optical_problem_stop(struct epitaxy *in);

//heat
long double epitaxy_get_heat_problem_start(struct epitaxy *in);
long double epitaxy_get_heat_problem_stop(struct epitaxy *in);
#endif
