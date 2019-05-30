// 
// General-purpose Photovoltaic Device Model gpvdm.com- a drift diffusion
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


/** @file data2d.h
	@brief Header file for data2d.c
*/
#ifndef data2d_h
#define data2d_h
#include "advmath.h"
#include <sim_struct.h>


struct data2d
{
	long double **data;
	int x_len;
	int y_len;
	long double *y_mesh;
	long double *x_mesh;

};

void data2d_init(struct data2d *in, int x_len, int y_len);
void data2d_set_value(struct data2d *in,long double value);
void data2d_free(struct data2d *in);
void data2d_init_y_mesh(struct data2d *in,long double start, long double stop);

#endif
