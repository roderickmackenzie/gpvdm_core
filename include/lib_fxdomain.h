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

/** @file lib_fxdomain.h
@brief Code to read input files.
*/

#ifndef lib_fxdomain_h
#define lib_fxdomain_h
#include "advmath.h"
#include "inp_struct.h"
#include <sim_struct.h>
#include "list_struct.h"

struct fit_sin_data
{
	int change;
	struct istruct orig_data;
	struct istruct guess_data;
	struct istruct delta_data;

	char *prefix;
	long double fx;

	double corse_mag;		//first guess at mag and delta
	double corse_delta;

	double delta;
	double mag;
	double error;
	struct simulation *sim;
};

struct fxdomain
{
	int fxdomain_sim_mode;
	int fxdomain_points;
	int fxdomain_n;
	long double fxdomain_Vexternal;
	long double fxdomain_voltage_modulation_max;
	long double fxdomain_light_modulation_depth;
	char fxdomain_modulation_type[100];
	int fxdomain_measure;
	long double fxdomain_L;

	//roll off
	int fxdomain_modulation_rolloff_enable;
	long double fxdomain_modulation_rolloff_start_fx;
	long double fxdomain_modulation_rolloff_speed;

	//dump verbosity
	int fxdomain_dump_verbocity;
	int fxdomain_screen_verbocity;
	char snapshot_path[PATH_MAX];
	char prefix_result[20];
	char prefix_modulation[20];
	long double fx;
	int modulate_voltage;
	int total_steps;

	//output data
	struct istruct out_j;
	struct istruct out_j_cut;

	struct istruct out_v;
	struct istruct out_v_cut;

	struct istruct out_modulation;
	struct istruct out_modulation_cut;

	//fitting
	int fxdomain_do_fit;
	int periods_to_fit;
	long double last_real;
	long double last_imag;
	long double last_mag;
	long double last_phi;
	long double last_j_error;
	long double last_v_error;
	long double last_mod_error;
	int large_signal;


};

//fitting
void fit_sin(struct simulation *sim,int dump_fit_progress_data,long double *fit_error,gdouble *ret_mag,gdouble *ret_delta,struct istruct *input_data,gdouble fx,char * prefix, char *output_path);
double sin_f (double *p,int len);
void fit_sin_dump(struct simulation *sim, struct fit_sin_data *fit_data, char *output_path);

//fxdomain
void fxdomain_solve(struct simulation *sim,struct device *in,struct fxdomain *fxdomain_config);
void fxdomain_init(struct simulation *sim,struct fxdomain *fxdomain_config);
void fxdomain_malloc(struct simulation *sim,struct fxdomain *fxdomain_config);
void fxdomain_dump(struct simulation *sim,struct fxdomain *fxdomain_config);
int fxdomain_fit(struct simulation *sim,struct device *in,struct fxdomain *fxdomain_config);
void fxdomain_free(struct simulation *sim,struct fxdomain *fxdomain_config);
void fxdomain_reset(struct simulation *sim,struct fxdomain *fxdomain_config);
void fxdomain_load_config(struct simulation *sim,struct fxdomain *fxdomain_config,struct device *dev,char *config_file_name);
void fxdomain_large_signal_solve(struct simulation *sim,struct device *in,struct fxdomain *fxdomain_config);
void fxdomain_small_signal_solve(struct simulation *sim,struct device *in,struct fxdomain *fxdomain_config);
void fxdomain_cal_complex_j(struct simulation *sim,struct device *in,struct newton_state_complex *ns,int z,int x);
#endif
