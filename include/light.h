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
#include <epitaxy_struct.h>
#include <ray.h>
#include <matrix.h>
#include <object.h>
#include <dat_file_struct.h>
#include <dim.h>



struct light
{
	char config_file[300];
	char dump_dir[PATH_MAX];
	struct dim_light dim;

	//zxyl
	long double ****Ep;
	long double ****Epz;
	long double ****En;
	long double ****Enz;
	long double ****n;
	long double ****alpha0;
	long double ****alpha;
	long double ****photons;
	long double ****photons_asb;
	long double ****pointing_vector;
	long double ****E_tot_r;
	long double ****E_tot_i;
	long double ****H;

	//complex zxyl	
	long double complex ****t;
	long double complex ****r;
	long double complex ****nbar;

	//3D arrrays
	long double ***Gn;
	long double ***Gp;
	long double ***Htot;
	long double ***photons_tot;
	struct object ****obj;

	//1D arrays
	long double *reflect;
	long double *transmit;

	//Input spectra
	struct istruct sun_read;
	long double *sun;
	long double *sun_norm;
	long double *sun_photons;
	long double *sun_E;
	char suns_spectrum_file[200];
	char light_file_generation[300];
	
	struct matrix mx;

	//laser
	long double laser_wavelength;
	int laser_pos;
	long double ND;
	long double spotx;
	long double spoty;
	long double pulseJ;
	long double pulse_width;

	long double *G_percent;
	//long double device_ylen;
	long double Psun;
	long double laser_eff;
	long double simplephotondensity;
	long double simple_alpha;
	long double Dphotoneff;

	//Dll section
	void (*fn_init)();
	void (*fn_solve_and_update)();
	int (*fn_solve_lam_slice)();
	long double (*fn_cal_photon_density)();
	void (*light_ver)();
	void *lib_handle;

	//config
	long double lstart;
	long double lstop;
	char mode[20];
	long double electron_eff;
	long double hole_eff;
	int force_update;

	long double *extract_eff;

	//Config values
	int align_mesh;
	int flip_field;
	int disable_transfer_to_electrical_mesh;
	int disable_cal_photon_density;
	long double light_file_generation_shift;
	long double light_flat_generation_rate;

	int print_wavlengths;
	int save_data_to_disk;

	int finished_solveing;

	long double last_Psun;
	long double last_laser_eff;
	long double last_wavelength_laser;


	struct epitaxy *epi;
};

#endif
