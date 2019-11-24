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

struct epi_layer
{
	int layer_number;
	long double y_start;
	long double y_stop;
	struct shape shapes[10];
	int nshape;
	char name[100];
	long double width;
	char pl_file[100];
	int pl_use_experimental_emission_spectra;
	long double pl_experimental_emission_efficiency;
	int pl_enabled;
	long double pl_fe_fh;
	long double pl_fe_te;
	long double pl_te_fh;
	long double pl_th_fe;
	long double pl_fh_th;
	char pl_spectrum_file[PATH_MAX];
	struct istruct pl_spectrum;
	double *photon_extract_eff;
	double *photon_extract_eff_count;
	long double avg_photon_extract_eff;
	long double peak_wavelength;

	long double shunt;
	long double series;
	long double C;
	long double n_ideality;
	long double J0;

	int electrical_layer;

	struct istruct alpha;
	struct istruct n;
};

struct epitaxy
{
	int layers;
	int electrical_layers;
	struct epi_layer layer[20];
	long double device_start;
	long double y_pos[20];
	char mat_file[20][100];
	char dos_file[20][100];

	char lumo_file[20][100];
	char homo_file[20][100];
	char shape_file[20][100];
	long double rgb[20][3];

};

void epitaxy_load(struct simulation *sim,struct epitaxy *in, char *file);
gdouble epitaxy_get_electrical_length(struct epitaxy *in);
gdouble epitaxy_get_optical_length(struct epitaxy *in);
int epitaxy_get_optical_material_layer(struct epitaxy *in,gdouble pos);
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
void epitaxy_load_electrical_file(struct simulation *sim,char *file_name, struct epi_layer *layer);
#endif
