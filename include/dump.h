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

/** @file dump.h
	@brief Fucntions to write data to disk.
*/


#ifndef h_dump
#define h_dump
#include "device.h"
#include "dump_ctrl.h"
#include "dynamic_store.h"
#include "dat_file.h"
#include <sim_struct.h>

void dump_init(struct simulation *sim,struct device* in);
void dump_load_config(struct simulation* sim,struct device *in);
void dump_remove_snapshots(struct simulation* sim);
void dump_dynamic_init(struct simulation *sim,struct dynamic_store *store,struct device *in);
void dump_dynamic_save(struct simulation *sim,struct device *in,char *outputpath,struct dynamic_store *store);
void dump_dynamic_add_data(struct simulation *sim,struct dynamic_store *store,struct device *in, gdouble x_value);
void dump_dynamic_free(struct simulation *sim,struct device *in,struct dynamic_store *store);
void dump_slice(struct device *in,char *prefix);
void dump_energy_slice(struct simulation *sim,char *out_dir,struct device *in);
void dump_device_map(struct simulation *sim,char* out_dir,struct device *in);
void dump_1d_slice(struct simulation *sim,struct device *in,char *out_dir);
void dump_write_to_disk(struct simulation *sim,struct device* in);
void buffer_add_3d_device_data(struct simulation *sim,struct dat_file *buf,struct device *in, gdouble ***data);
void buffer_set_graph_type(struct dat_file *buf,struct device *in);
void dumpfiles_load(struct simulation* sim);
void dumpfiles_free(struct simulation* sim);
int dumpfiles_should_dump(struct simulation* sim,char *name);
void dumpfiles_process(struct simulation* sim,struct istruct *in,char *name);
void dumpfiles_turn_on_dir(struct simulation* sim,char *in);
void dump_clean_cache_files(struct simulation* sim);
void dump_make_snapshot_dir(struct simulation *sim,char *out_dir ,long double time,long double voltage, long double fx, int number);
void dump_make_snapshot_dir_with_name(struct simulation *sim,char *out_dir ,long double time,long  double voltage, long double fx, int number,char *snapshot_name);
void buffer_add_3d_device_data_int(struct simulation *sim,struct dat_file *buf,struct device *in,int ***data);
void buffer_add_3d_to_2d_projection(struct simulation *sim,struct dat_file *buf,struct device *in,gdouble ***data);
#endif
