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


/** @file perovskite.h
	@brief Header file for perovskite.c
*/
#ifndef perovskite_h
#define perovskite_h

#include <sim_struct.h>
#include <device.h>
#include <perovskite_struct.h>

int perovskite_enabled(struct device *in);
void perovskite_dump(struct simulation *sim,struct device *in);
long double get_perovskite_error(struct device *in,long double *b);
int solve_perovskite(struct simulation *sim,struct device *in, int z, int x);
void perovskite_build_matrix(struct simulation *sim,struct device *in,struct perovskite *config, int z, int x);
void perovskite_init(struct simulation *sim,struct device *in);
void perovskite_free(struct simulation *sim,struct device *in);


#endif
