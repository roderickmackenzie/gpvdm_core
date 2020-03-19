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

/** @file device.h
	@brief The main structure which holds information about the device.
*/

#ifndef device_h
#define device_h
#include <stdio.h>
#include "code_ctrl.h"
#include "light.h"
#include <epitaxy_struct.h>
#include "advmath.h"
#include <dos_struct.h>
#include <contact_struct.h>
#include <perovskite_struct.h>
#include <circuit_struct.h>
#include <dim.h>
#include <matrix.h>
#include <shape_struct.h>
#include <heat.h>
#include <mesh_struct.h>


struct solver_cache
{
	char hash[100];
	int length;
	int enabled;
};


struct newton_state
{
	struct dimensions dim;

	int last_ittr;
	long double last_error;

	long double ***phi;
	long double ***x;
	long double ***xp;

	long double ****xt;
	long double ****xpt;

};

struct newton_state_complex
{
	struct dimensions dim;

	long double complex ***phi;
	long double complex ***x;
	long double complex ***xp;

	long double complex ****xt;
	long double complex ****xpt;

};

struct device
{
	struct epitaxy my_epitaxy;
	//Device state
		//0D arrays

		//mesh points

		int remesh;
		int newmeshsize;

		int dynamic_mesh;

		int excite_conv;

		gdouble deltaFln;
		gdouble deltaFlp;
		gdouble deltaFrn;
		gdouble deltaFrp;

		gdouble xlen;
		gdouble ylen;
		gdouble zlen;

		//1D arrays

		//2D arrays
		long double **Vapplied_y0;
		long double **Vapplied_y1;
		long double **Vapplied_x0;
		long double **Vapplied_x1;

		int **passivate_y0;
		int **passivate_y1;
		int **passivate_x0;
		int **passivate_x1;

		int **n_contact_y0;
		int **n_contact_y1;
		int **n_contact_x0;
		int **n_contact_x1;

		long double **Jn_y0;
		long double **Jn_y1;
		long double **Jp_y0;
		long double **Jp_y1;

		long double **Fi0_y0;		//This is the equilibrium fermi level of the contact were it in free space, i.e. with no phi subtracted
		long double **Fi0_y1;		//This is referenced to Fi0_y0[0][0], and is the difference between Fi0_y0[0][0] and Fi0_y0[z][x/y], this difference must be equal to the built in potential on the contact.
		long double **Fi0_x0;
		long double **Fi0_x1;

		long double **V_y0;
		long double **V_y1;
		long double **V_x0;
		long double **V_x1;

		//3D arrays zxy
		gdouble ***n;
		gdouble ***p;
		gdouble ***dn;
		gdouble ***dndphi;
		gdouble ***dp;
		gdouble ***dpdphi;


		gdouble ***Rfree;

		//material constants
		gdouble ***mun;
		gdouble ***mup;
		gdouble ***Nad;
		gdouble ***muion;
		gdouble ***Eg;
		gdouble ***Xi;
		gdouble ***Ev;
		gdouble ***Ec;

		gdouble ***G;
		gdouble ***Gn;
		gdouble ***Gp;
		gdouble ***Photon_gen;

		gdouble ***Dn;
		gdouble ***Dp;
		gdouble ***epsilonr;

		int ***imat;

		//Calculated

		gdouble ***Nion;
		gdouble ***Nion_last;


		gdouble ***Fn;
		gdouble ***Fp;
		gdouble ***Nc;
		gdouble ***Nv;

		gdouble ***Fi;

		int ***imat_epitaxy;
		int ***mask;

		gdouble ***Jn;
		gdouble ***Jp;
		gdouble ***Jn_x;
		gdouble ***Jp_x;

		gdouble ***Jn_diffusion;
		gdouble ***Jn_drift;

		gdouble ***Jp_diffusion;
		gdouble ***Jp_drift;


		//gdouble ***x;
		gdouble ***t;
		//gdouble ***xp;
		gdouble ***tp;
		gdouble ***kf;
		gdouble ***kd;
		gdouble ***kr;

		gdouble ***Rn;
		gdouble ***Rp;
		gdouble ***Rn_srh;
		gdouble ***Rp_srh;
		gdouble ***Rnet;

		gdouble ***Rbi_k;

		gdouble ***ex;
		gdouble ***Dex;
		gdouble ***Hex;

		gdouble ***nf_save;
		gdouble ***pf_save;
		gdouble ***nt_save;
		gdouble ***pt_save;

		gdouble ***nfequlib;
		gdouble ***pfequlib;
		gdouble ***ntequlib;
		gdouble ***ptequlib;

		gdouble ***phi_save;

		gdouble ***nlast;
		gdouble ***plast;

		gdouble ***wn;
		gdouble ***wp;

		//n traps
			gdouble ***nt_all;
			gdouble ***tt;

		//p traps
			gdouble ***pt_all;
			gdouble ***tpt;


		gdouble ***nrelax;
		gdouble ***ntrap_to_p;
		gdouble ***prelax;
		gdouble ***ptrap_to_n;

		gdouble ***n_orig;
		gdouble ***p_orig;
		gdouble ***n_orig_f;
		gdouble ***p_orig_f;
		gdouble ***n_orig_t;
		gdouble ***p_orig_t;

		gdouble ***B;

		//4D arrays
		gdouble ****ntb_save;
		gdouble ****ptb_save;


		//n traps
			gdouble ****nt;
			gdouble ****ntlast;
			gdouble ****dnt;
			gdouble ****srh_n_r1;
			gdouble ****srh_n_r2;
			gdouble ****srh_n_r3;
			gdouble ****srh_n_r4;
			gdouble ****dsrh_n_r1;
			gdouble ****dsrh_n_r2;
			gdouble ****dsrh_n_r3;
			gdouble ****dsrh_n_r4;
			gdouble ****Fnt;
			//gdouble ****xt;


			gdouble ****nt_r1;
			gdouble ****nt_r2;
			gdouble ****nt_r3;
			gdouble ****nt_r4;
		//p traps
			gdouble ****pt;
			gdouble ****ptlast;
			gdouble ****dpt;
			gdouble ****srh_p_r1;
			gdouble ****srh_p_r2;
			gdouble ****srh_p_r3;
			gdouble ****srh_p_r4;
			gdouble ****dsrh_p_r1;
			gdouble ****dsrh_p_r2;
			gdouble ****dsrh_p_r3;
			gdouble ****dsrh_p_r4;
			gdouble ****Fpt;
			//gdouble ****xpt;


			gdouble ****pt_r1;
			gdouble ****pt_r2;
			gdouble ****pt_r3;
			gdouble ****pt_r4;

	struct matrix mx;

	//Arrays used by newton solver
	gdouble *newton_dntrap;
	gdouble *newton_dntrapdntrap;
	gdouble *newton_dntrapdn;
	gdouble *newton_dntrapdp;
	gdouble *newton_dJdtrapn;
	gdouble *newton_dJpdtrapn;

	gdouble *newton_dptrapdp;
	gdouble *newton_dptrapdptrap;
	gdouble *newton_dptrap;
	gdouble *newton_dptrapdn;
	gdouble *newton_dJpdtrapp;
	gdouble *newton_dJdtrapp;
	gdouble *newton_dphidntrap;
	gdouble *newton_dphidptrap;
	gdouble *newton_ntlast;
	gdouble *newton_ptlast;

	//math
	int max_electrical_itt;
	gdouble electrical_clamp;
	int max_electrical_itt0;
	gdouble electrical_clamp0;
	gdouble electrical_error0;
	int math_enable_pos_solver;
	gdouble min_cur_error;
	int pos_max_ittr;
	char solver_name[20];
	char complex_solver_name[20];
	char newton_name[20];

	//meshing
	struct mesh_obj mesh_data;
	struct mesh_obj mesh_data_save;


	gdouble dt;
	int srh_sim;
	int go_time;
	gdouble time;
	long double fx;


	int stop;
	gdouble Rshort;



	int onlypos;
	int odes;
	gdouble posclamp;

	gdouble A;
	gdouble Vol;

	gdouble Rshunt;
	gdouble Rcontact;
	gdouble Rload;
	gdouble L;
	gdouble C;

	int stop_start;
	gdouble externalv;
	gdouble Ilast;
	int timedumpcount;
	char simmode[200];
	gdouble area;

	int ntrapnewton;
	int ptrapnewton;

	long double **electrons_y0;
	long double **holes_y0;

	long double **electrons_y1;
	long double **holes_y1;

	long double **electrons_x0;
	long double **holes_x0;

	long double **electrons_x1;
	long double **holes_x1;


	gdouble t_big_offset;

	gdouble other_layers;

	int kl_in_newton;
	int config_kl_in_newton;
	void (*newton_aux)(struct simulation *sim, struct device* ,gdouble ,gdouble* ,gdouble* ,gdouble* ,gdouble* ,gdouble* ,gdouble* ,gdouble* ,gdouble*);
	gdouble xnl_left;
	gdouble xpl_left;
	int stoppoint;
	gdouble ilast;

	int newton_clever_exit;
	char plot_file[100];

	gdouble start_stop_time;


	gdouble Is;
	gdouble n_id;
	gdouble Igen;

	//Light
	struct light mylight;
	struct light probe_modes;
	struct math_xy steady_stark;

	gdouble Vbi;
	int newton_min_itt;
	gdouble vbi;
	gdouble avg_gen;

	//dump
	int dump_energy_slice_xpos;
	int dump_energy_slice_ypos;
	int dump_energy_slice_zpos;

	int dump_1d_slice_xpos;
	int dump_1d_slice_zpos;

	int dumpitdos;

	long double dump_dynamic_pl_energy;

	int snapshot_number;

	//pl
	long double pl_intensity;
	long double pl_intensity_tot;
	int pl_enabled;
	int pl_use_experimental_emission_spectra;
	gdouble Rext;
	gdouble Cext;
	gdouble VCext_last;
	gdouble VCext;
	int newton_last_ittr;
	gdouble phi_mul;
	long double layer_start[100];
	long double layer_stop[100];


	struct dos dosn[10];
	struct dos dosp[10];

	gdouble *tm_sun;
	gdouble *tm_voltage;
	gdouble *tm_laser;
	gdouble *tm_time_mesh;
	gdouble *tm_fs_laser;
	int tm_mesh_len;
	int tm_use_mesh;
	int tm_mesh_pos;
	int dd_conv;

	//thermal
	struct heat thermal;
	long double ***Tl;
	long double ***Te;
	long double ***Th;

	long double ***Hl;
	long double ***H_recombination;
	long double ***H_joule;

	long double ***He;
	long double ***Hh;

	long double ***ke;
	long double ***kh;

	//contacts
	struct contact contacts[10];
	int ncontacts;
	int active_contact;

	struct perovskite mobileion;

	gdouble map_start;
	gdouble map_stop;

	long double flip_current;

	//Ray tracing

	struct image my_image;

	//solver cache
	struct solver_cache cache;
	int newton_only_fill_matrix;

	struct newton_state ns;
	struct newton_state ns_save;

	long double omega;

	struct circuit cir;
	int circuit_simulation;

	//objects
	struct object *obj;		//This is the scene built from triangles
	int objects;
	int triangles;

	struct shape big_box;

};


#endif
