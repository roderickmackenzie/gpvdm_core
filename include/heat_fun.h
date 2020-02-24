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

/** @file heat_fun.h
@brief heat functions from libheat
*/

#ifndef h_heat_fun
#define h_heat_fun
#include <complex.h>
#include "advmath.h"
#include "i.h"
#include <sim_struct.h>
#include <epitaxy_struct.h>
#include <heat.h>
#include <device.h>

//General
int heat_solve(struct simulation *sim, struct heat *thermal,struct device *dev, int z, int x);
void heat_load_config(struct simulation *sim,struct heat *thermal, struct device *dev);
void heat_load_config_file(struct simulation *sim,struct heat *thermal);
void heat_setup_dump_dir(struct simulation *sim,struct heat *thermal);
void heat_dump(struct simulation *sim,struct heat *thermal);
void heat_free_memory(struct simulation *sim,struct heat *thermal);
void heat_malloc(struct simulation *sim,struct heat *thermal);
void heat_init(struct heat *thermal);

//Transfer between meshes
void heat_transfer_device_heat_to_heat_mesh(struct heat *thermal, struct device *dev);
void heat_transfer_optical_heat_to_heat_mesh(struct heat *thermal, struct light *li);

//New thermal model
void heat_build_obj_pointer_array(struct simulation *sim,struct heat *thermal, struct device *dev);
void heat_build_materials_arrays(struct simulation *sim,struct heat *thermal, struct device *dev);

//Hydrodynamic model
void hydrodynamic_transfer_temperatures_to_device(struct simulation *sim,struct device *dev,struct heat *thermal,int z, int x);
void hydrodynamic_update_heat(struct simulation *sim, struct heat *thermal,struct device *dev,int z, int x);
void hydrodynamic_mesh_check(struct simulation *sim, struct heat *thermal,struct device *dev);
int hydrodynamic_solve(struct simulation *sim, struct heat *thermal,struct device *dev, int z, int x);


#endif
