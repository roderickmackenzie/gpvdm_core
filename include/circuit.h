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

/** @file circuit.h
@brief Header files for nodal analysis
*/

#ifndef circuit_h
#define circuit_h
#include <sim_struct.h>
#include <circuit_struct.h>

void circuit_init(struct circuit *cir);
void circuit_malloc(struct simulation * sim,struct circuit *cir);
void circuit_free(struct simulation * sim,struct circuit *cir);
void circuit_load(struct simulation * sim,struct circuit *cir);
void circuit_alloc_nodes_and_links(struct simulation * sim,struct circuit *cir);
void circuit_build_device(struct simulation * sim,struct circuit *cir,struct device *dev);
void circuit_solve(struct simulation * sim,struct circuit *cir,struct device *dev);
void circuit_apply_voltages(struct simulation * sim,struct circuit *cir,struct device *dev);
void circuit_transfer_to_electrical_mesh(struct simulation * sim,struct circuit *cir,struct device *dev);
double circuit_node_get_j(struct simulation * sim,struct circuit *cir,int n);
void circuit_printf_links(struct simulation * sim,struct circuit *cir);
void circuit_dump_j(struct simulation * sim,struct device *dev,char *out_dir);
struct circuit_node *circuit_add_node(struct simulation * sim,struct circuit *cir,int x,int y,int z,int type);
int circuit_find_node_by_xyz(struct simulation * sim,struct circuit *cir,int x,int y,int z);
void circuit_node_set_type(struct simulation * sim, struct circuit_node *node, int type);
void circuit_print_nodes(struct simulation * sim,struct circuit *cir);
void circuit_print_links(struct simulation * sim,struct circuit *cir);
struct circuit_link *circuit_add_link(struct simulation * sim,struct circuit *cir,struct circuit_link *in_link);
int circuit_load_config(struct simulation * sim,struct circuit *cir);
void circuit_time_step(struct simulation * sim,struct circuit *cir);
void circuit_apply_photo_generation(struct simulation * sim,struct circuit *cir,struct device *dev);
void circuit_dump(struct simulation * sim,struct device *dev,struct circuit *cir);
void circuit_cal_resistance(struct simulation * sim,struct circuit *cir,struct device *dev);
void circuit_plot_resistance(struct simulation * sim,struct circuit *cir,struct device *dev);
int circuit_find_link(struct simulation * sim,struct circuit *cir,int z0,int x0,int y0, int z1,int x1,int y1);
double circuit_get_max_j(struct simulation * sim,struct circuit *cir);
#endif
