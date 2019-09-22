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

/** @file light.h
@brief light functions from liblight
*/

#ifndef h_light
#define h_light
#include <complex.h>
#include "advmath.h"
#include "i.h"
#include <sim_struct.h>
#include <epitaxy.h>
#include <ray.h>

struct light
{
	char config_file[300];
	char dump_dir[1024];
	int points;
	int lpoints;
	gdouble *x;
	gdouble dx;
	gdouble **Ep;
	gdouble **Epz;
	gdouble **En;
	gdouble **Enz;
	gdouble **n;
	gdouble **alpha0;
	gdouble **alpha;
	gdouble **photons;
	gdouble **photons_asb;
	gdouble **pointing_vector;
	gdouble *photons_tot;
	gdouble **E_tot_r;
	gdouble **E_tot_i;
	gdouble **H;
	gdouble *reflect;
	gdouble *transmit;
	int *layer;
	gdouble *sun_E;
	gdouble *H1d;
	gdouble ylen;

	struct istruct sun_read;
	gdouble *sun;
	gdouble *sun_norm;
	gdouble *sun_photons;
	char suns_spectrum_file[200];
	char light_file_generation[300];
	int M;
	int N;
	int *Ti;
	int *Tj;
	double *Tx;
	double *Txz;
	double *b;
	double *bz;
	gdouble lstart;
	gdouble lstop;
	gdouble *l;
	gdouble *Gn;
	gdouble *Gp;
	gdouble dl;
	gdouble laser_wavelength;
	int laser_pos;
	gdouble ND;
	gdouble spotx;
	gdouble spoty;
	gdouble pulseJ;
	gdouble pulse_width;
	gdouble complex **t;
	gdouble complex **r;
	gdouble complex **nbar;
	gdouble *layer_end;
	gdouble device_start;
	gdouble *G_percent;
	gdouble device_ylen;
	gdouble Eg;
	gdouble Psun;
	gdouble laser_eff;
	gdouble simplephotondensity;
	gdouble simple_alpha;
	gdouble Dphotoneff;
	void (*fn_init)();
	void (*fn_solve_and_update)();
	int (*fn_solve_lam_slice)();
	gdouble (*fn_cal_photon_density)();
	void (*light_ver)();
	void *lib_handle;
	char mode[20];
	gdouble electron_eff;
	gdouble hole_eff;
	int force_update;

	int device_start_layer;
	int device_start_i;

	long double *extract_eff;

	//Flags
	int align_mesh;
	int flip_field;
	int disable_transfer_to_electrical_mesh;
	int disable_cal_photon_density;
	long double light_file_generation_shift;

	//should remove this
	struct epitaxy *my_epitaxy;
	//
	int print_wavlengths;
	int save_data_to_disk;

	int finished_solveing;
};

#endif
