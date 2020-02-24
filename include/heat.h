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

/** @file heat.h
@brief a structure for the heat model
*/

#ifndef h_heat
#define h_heat
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
#include <mesh_struct.h>


struct heat
{
	char dump_dir[PATH_MAX];
	struct dim_heat dim;
	int thermal_model_type;

	//zxy
	long double ***Tl;
	long double ***Te;
	long double ***Th;

	long double ***Hl;
	long double ***He;
	long double ***Hh;

	long double ***kl;
	long double ***ke;
	long double ***kh;


	long double thermal_kl;

	long double thermal_tau_e;
	long double thermal_tau_h;

	struct matrix mx;

	struct object ****obj;

	//boundry conditions
	long double Ty0;
	long double Ty1;
	long double Tx0;
	long double Tx1;
	long double Tz0;
	long double Tz1;

	int Tliso;
	int Triso;
	int nofluxl;

	//convergance
	int thermal_conv;
	long double min_error;
	int newton_enable_external_thermal;
	int thermal_l;
	int thermal_e;
	int thermal_h;

	struct mesh_obj mesh_data;
};

#endif
