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

/** @file dynamic_store.h
	@brief Store information as the simulation progresses such as voltage, current or carrier density, these are 1D arrays as a function of time of simulation step.
*/

#ifndef dynamic_store_h
#define dynamic_store_h
#include "i.h"

struct dynamic_store
{
	//recombination
	struct istruct R_nfree_to_pfree;
	struct istruct R_srh_nfree;
	struct istruct R_srh_pfree;
	struct istruct R_srh_nfree_to_ptrap;
	struct istruct R_srh_pfree_to_ntrap;
	struct istruct T_srh_pfree_to_ptrap;
	struct istruct T_srh_nfree_to_ntrap;
	struct istruct G_n;
	struct istruct G_p;
	struct istruct R_surface_y0;
	struct istruct R_surface_y1;

	//charge
	struct istruct Q_nfree;
	struct istruct Q_pfree;
	struct istruct Q_ntrap;
	struct istruct Q_ptrap;
	struct istruct Q_nfree_and_ntrap;
	struct istruct Q_pfree_and_ptrap;

	//struct istruct charge_change;
	//struct istruct ntrap_delta_out;
	//struct istruct ptrap_delta_out;
	//struct istruct nfree_delta_out;
	//struct istruct pfree_delta_out;
	//struct istruct dynamic_charge_tot;
	//struct istruct tpc_filledn;
	//struct istruct tpc_filledp;
	//struct istruct dynamic_np;


	//mobility
	struct istruct mu_n;
	struct istruct mu_p;
	struct istruct mu_n_p_avg;

	//srh rates
	struct istruct srh_n_r1;
	struct istruct srh_n_r2;
	struct istruct srh_n_r3;
	struct istruct srh_n_r4;

	struct istruct srh_p_r1;
	struct istruct srh_p_r2;
	struct istruct srh_p_r3;
	struct istruct srh_p_r4;


	//J
	struct istruct J_y0_n;
	struct istruct J_y0_p;
	struct istruct J_y1_n;
	struct istruct J_y1_p;

	struct istruct jout;
	struct istruct jn_avg;
	struct istruct jp_avg;
	struct istruct dynamic_jn;
	struct istruct dynamic_jp;
	struct istruct jnout_mid;
	struct istruct jpout_mid;
	struct istruct iout;
	struct istruct dynamic_jn_drift;
	struct istruct dynamic_jn_diffusion;

	struct istruct dynamic_jp_drift;
	struct istruct dynamic_jp_diffusion;

	//pl
	struct istruct dynamic_pl;
	struct istruct dynamic_pl_tot;

	//field
	struct istruct E_field;
	struct istruct dynamic_Vapplied;
	struct istruct band_bend;

	//other
	struct istruct dynamic_qe;
	long double ***band_snapshot;
};
#endif
