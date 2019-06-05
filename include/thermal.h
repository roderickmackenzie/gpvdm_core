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

/** @file thermal.h
@brief header file for the theral solver, not really used yet.
*/

#ifndef thermal_h
#define thermal_h
#include "sim.h"
void update_heat(struct device *in);
void dump_thermal(struct simulation *sim,struct device *in);
int solve_thermal(struct simulation *sim,struct device *in);
void thermal_init(struct device *in);
void thermal_free(struct device *in);
double get_thermal_error(struct device *in);
double get_last_thermal_error();
#endif
