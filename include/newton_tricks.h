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

/** @file newton_tricks.h
@brief functions to solver for an external voltage
*/

#ifndef newton_tricks_h
#define newton_tricks_h

struct newton_math_state
{
	int max_electrical_itt;
	gdouble min_cur_error;
	int newton_min_itt;
	gdouble electrical_clamp;
	int newton_clever_exit;
};


void newton_push_state(struct device *in);
void newton_pop_state(struct device *in);
gdouble sim_externalv(struct simulation *sim,struct device *in,gdouble wantedv);
gdouble sim_i(struct simulation *sim,struct device *in,gdouble wantedi);
void auto_ramp_contacts(struct simulation *sim,struct device *in);
void ramp_externalv(struct simulation *sim,struct device *in,gdouble from,gdouble to);
void set_ntricks_fast(int val);
gdouble sim_voc(struct device *in);
void newton_sim_simple(struct simulation  *sim,struct device *in);
void ntricks_auto_ramp_contacts(struct simulation *sim,struct device *in);


void newton_externv_aux(struct simulation *sim,struct device *in,gdouble V,gdouble* i,gdouble* didv,gdouble* didphi,gdouble* didxil,gdouble* didxipl,gdouble* didphir,gdouble* didxir,gdouble* didxipr);
gdouble newton_externv(struct simulation *sim,struct device *in,gdouble Vtot,int usecap);
gdouble newton_externalv_simple(struct simulation *sim,struct device *in,gdouble V);
long double sim_externalv_ittr(struct simulation *sim,struct device *in,gdouble wantedv);

//perovskite functions
long double newton_externalv_simple_perovskite(struct simulation *sim,struct device *in,gdouble V);
long double newton_externv_perovskite(struct simulation *sim,struct device *in,gdouble Vtot,int usecap);

void state_cache_init(struct simulation *sim,struct device *in);
void hash_dir(struct simulation *sim,char *out);
void state_gen_vector(struct simulation *sim,struct device *in);
int state_find_vector(struct simulation *sim,struct device *in,char *out);

//newton state
void newton_state_init(struct newton_state *ns);
void newton_state_alloc_mesh(struct newton_state *ns,struct dimensions *dim);
void newton_state_alloc_traps(struct newton_state *ns,struct dimensions *dim);
void newton_state_free(struct newton_state *ns);
void newton_state_cpy(struct newton_state *out,struct newton_state *in);
void newton_state_save(struct simulation *sim,char *file_name,struct newton_state *ns);
int newton_state_load(struct simulation *sim,struct newton_state *ns,char *file_name);
void newton_state_update_device(struct simulation *sim,struct device *in, struct newton_state *ns);

//newton state complex
void newton_state_complex_init(struct newton_state_complex *ns);
void newton_state_complex_alloc_mesh(struct newton_state_complex *ns,struct dimensions *dim);
void newton_state_complex_alloc_traps(struct newton_state_complex *ns,struct dimensions *dim);
void newton_state_complex_free(struct newton_state_complex *ns);
#endif

