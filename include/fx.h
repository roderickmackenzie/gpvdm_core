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
/** @file fx.h
@brief fxdomain solver
*/


#ifndef _fx
#define _fx
#include <sim_struct.h>

void fx_mesh_save(struct simulation *sim);
void fx_load_mesh(struct simulation *sim,struct device *in,int number);
void fx_step(struct simulation *sim,struct device *in);
void fxdomain_load_config(struct simulation *sim,struct fxdomain *in,struct device *dev,char *config_file_name);
int fx_run();
long double fx_get_fx();
void fx_memory_free();
int fx_points();
int fx_get_step();
#endif
