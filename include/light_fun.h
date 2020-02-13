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

#ifndef h_light_fun
#define h_light_fun
#include <complex.h>
#include "advmath.h"
#include "i.h"
#include <sim_struct.h>
#include <epitaxy_struct.h>
#include <light.h>
#include <device.h>

void light_norm_photon_density(struct simulation *sim, struct light *li);
void light_malloc(struct simulation *sim,struct light *li);
void light_load_materials(struct simulation *sim,struct light *li);
long double light_cal_photon_density(struct simulation *sim,struct light *li);
void light_load_config(struct simulation *sim,struct light *li, struct device *dev);
void light_load_config_file(struct simulation *sim,struct light *li);
void light_build_mesh(struct simulation *sim,struct light *li,struct dimensions *electrical_dim);
void light_build_materials_arrays(struct simulation *sim,struct light *li, struct device *dev);
void light_set_sun_power(struct light *li,long double power, long double laser_eff);
void light_free_memory(struct simulation *sim,struct light *li);
void light_get_mode(struct istruct *mode,int lam,struct light *li);
void light_set_unity_power(struct light *li);
void light_solve_optical_problem(struct simulation *sim,struct device *dev,struct light *li);
void light_solve_all(struct simulation *sim,struct device *dev,struct light *li);
void light_set_dump(struct light *li,int dump);
void light_free(struct simulation *sim,struct light *li);
void light_dump(struct simulation *sim,struct light *li);
int light_solve_lam_slice(struct simulation *sim,struct device *dev, struct light *li,int z, int x,int lam);
void light_dump_1d(struct simulation *sim,struct light *li, int i,char *ext);
void light_dump_verbose_1d(struct simulation *sim,struct light *li, int i,char *ext);
void light_dump_verbose_2d(struct simulation *sim,struct light *li);
void light_get_mode(struct istruct *mode,int lam,struct light *li);
void light_set_unity_laser_power(struct light *li,int lam);
void light_free_epitaxy(struct light *li);
void light_import_epitaxy(struct simulation *sim,struct light *li,struct epitaxy *my_epitaxy);
void light_calculate_complex_n(struct light *li);
int light_load_laser(struct simulation *sim, struct light *li,char *name);
long double light_get_sun(struct light *li);
void light_set_sun(struct light *li,long double Psun);
void light_set_model(struct light *li,char *model);
void light_set_sun_delta_at_wavelength(struct simulation *sim,struct light *li,long double lam);
void light_free_dlls(struct simulation *sim,struct light *li);
int light_get_pos_from_wavelength(struct simulation *sim,struct light *li,double lam);
void light_setup_dump_dir(struct simulation *sim,struct light *li);
long double light_calculate_photons_absorbed_in_active_layer(struct light *li);
void light_dump_sim_info(struct simulation *sim,struct light *li);
void light_transfer_gen_rate_to_device(struct device *dev,struct light *li);
#endif
