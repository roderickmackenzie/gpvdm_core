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

/** @file newton_interface.h
@brief an interface to the newton solvers.
*/


#ifndef h_newton_interface
#define h_newton_interface
#include <sim_struct.h>

void newton_set_min_ittr(struct device *in,int ittr);
void newton_init(struct simulation *sim,char *solver_name);
int solve_cur(struct simulation *sim,struct device *in,int z, int x);
void solver_realloc(struct simulation *sim,struct device * in);
void solver_free_memory(struct simulation *sim,struct device * in);
void newton_interface_free(struct simulation *sim);
#endif
