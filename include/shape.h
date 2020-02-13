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


/** @file epitaxy.h
	@brief Read the epitaxy from the epitaxy.inp file.
*/


#ifndef shape_h
#define shape_h
#include "advmath.h"
#include <sim_struct.h>
#include <shape_struct.h>

struct shape *shape_load_file(struct simulation *sim,struct epitaxy *in,struct shape *s, char *file_name);
int shape_get_index(struct simulation *sim,struct epitaxy *in,long double x,long double y,long double z);
void shape_free_materials(struct simulation *sim,struct epitaxy *in);
void shape_free(struct simulation *sim,struct shape *s);
int shape_in_shape(struct simulation *sim,struct shape *s,long double z,long double x,long double y);
void shape_init(struct simulation *sim,struct shape *s);
void shape_load_materials(struct simulation *sim,struct shape *s);
#endif
