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

#ifndef device_fun_h
#define device_fun_h
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
#include <device.h>
#include <shape.h>

void device_init(struct simulation *sim,struct device *in);
void device_alloc_traps(struct device *in);
void device_get_memory(struct simulation *sim,struct device *in);
void device_free(struct simulation *sim,struct device *in);
void device_load_math_config(struct simulation *sim,struct device *in);
void device_dump_world_to_file(struct simulation *sim,struct device *dev,char *file_name);
void device_build_scene(struct simulation *sim,struct device *dev);
void device_add_shape_to_world(struct simulation *sim,struct device *dev,struct shape *s);
void device_calculate_joule_heat(struct simulation *sim,struct device *dev);
void device_calculate_recombination_heat(struct simulation *sim,struct device *dev);

#endif
