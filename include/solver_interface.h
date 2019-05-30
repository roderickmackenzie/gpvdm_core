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

/** @file solver_interface.h
@brief an interface load and run sparse solvers from various dlls.
*/

#ifndef h_solver_interface
#define h_solver_interface
#include <sim_struct.h>
#include <device.h>
void solver_init(struct simulation *sim,char *solver_name);
void solver(struct simulation *sim,int col,int nz,int *Ti,int *Tj, long double *Tx,long double *b);
void dump_matrix(struct simulation *sim,struct device *in);
void solver_free(struct simulation *sim);

#endif
