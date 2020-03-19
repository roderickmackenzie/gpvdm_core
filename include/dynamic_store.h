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
	struct math_xy R_nfree_to_pfree;
	struct math_xy R_srh_nfree;
	struct math_xy R_srh_pfree;
	struct math_xy R_srh_nfree_to_ptrap;
	struct math_xy R_srh_pfree_to_ntrap;
	struct math_xy T_srh_pfree_to_ptrap;
	struct math_xy T_srh_nfree_to_ntrap;
	struct math_xy G_n;
	struct math_xy G_p;
	struct math_xy R_surface_y0;
	struct math_xy R_surface_y1;

	//charge
	struct math_xy Q_nfree;
	struct math_xy Q_pfree;
	struct math_xy Q_ntrap;
	struct math_xy Q_ptrap;
	struct math_xy Q_nfree_and_ntrap;
	struct math_xy Q_pfree_and_ptrap;

	//struct math_xy charge_change;
	//struct math_xy ntrap_delta_out;
	//struct math_xy ptrap_delta_out;
	//struct math_xy nfree_delta_out;
	//struct math_xy pfree_delta_out;
	//struct math_xy dynamic_charge_tot;
	//struct math_xy tpc_filledn;
	//struct math_xy tpc_filledp;
	//struct math_xy dynamic_np;


	//mobility
	struct math_xy mu_n;
	struct math_xy mu_p;
	struct math_xy mu_n_p_avg;

	//srh rates
	struct math_xy srh_n_r1;
	struct math_xy srh_n_r2;
	struct math_xy srh_n_r3;
	struct math_xy srh_n_r4;

	struct math_xy srh_p_r1;
	struct math_xy srh_p_r2;
	struct math_xy srh_p_r3;
	struct math_xy srh_p_r4;


	//J
	struct math_xy J_y0_n;
	struct math_xy J_y0_p;
	struct math_xy J_y1_n;
	struct math_xy J_y1_p;

	struct math_xy jout;
	struct math_xy jn_avg;
	struct math_xy jp_avg;
	struct math_xy dynamic_jn;
	struct math_xy dynamic_jp;
	struct math_xy jnout_mid;
	struct math_xy jpout_mid;
	struct math_xy iout;
	struct math_xy dynamic_jn_drift;
	struct math_xy dynamic_jn_diffusion;

	struct math_xy dynamic_jp_drift;
	struct math_xy dynamic_jp_diffusion;

	//pl
	struct math_xy dynamic_pl;
	struct math_xy dynamic_pl_tot;

	//field
	struct math_xy E_field;
	struct math_xy dynamic_Vapplied;
	struct math_xy band_bend;

	//other
	struct math_xy dynamic_qe;
	long double ***band_snapshot;
};
#endif
