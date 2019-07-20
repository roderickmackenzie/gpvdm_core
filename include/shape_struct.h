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


/** @file shape_struct.h
	@brief A structure to hold shapes
*/


#ifndef shape_struct_h
#define shape_struct_h
#include "advmath.h"
#include <sim_struct.h>

struct shape
{
	long double dx;
	long double dy;
	long double dz;
	long double dx_padding;
	long double dy_padding;
	long double dz_padding;
	int nx;
	int nz;
	char name[20];
	char dos_file[20];
	char shape_type[20];
	char optical_material[100];
	long double x0;
	long double y0;
	long double z0;
	int dos_index;
	int epi_index;
	struct istruct alpha;
	struct istruct n;
};

#endif
