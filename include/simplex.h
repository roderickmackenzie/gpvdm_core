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

/** @file simplex.c
@brief Downhill simplex code
*/
#ifndef simplex_h
#define simplex_h

#define TINY 1.0e-10		//A small number.
#define NMAX 5000			//Maximum allowed number of function evaluations.
#define SIMPLEX_CONVERGED 1
#define  SIMPLEX_MAX 2

#include <sim_struct.h>

struct multimin
{
	int ittr;
	int n_max;
	int nsimplex;
	int ndim;
	double stop_error;
	double *x;
	double **p;
	int i_hi0;
	int i_hi1;
	int i_lo;
	double *y;
	double *center;
	double  *ptry;
	double ytry;
	double *s;
	double error;
	double error_delta;
	double error_last;
	double (*fn)(double *p,int len);
};

void multimin_dump(struct simulation *sim,struct multimin *data);
void multimin_init(struct multimin *data);
void multimin_malloc(struct multimin *data);
void multimin_init_simplex(struct multimin *data);
void multimin_free(struct multimin *data);
void multimin_cal_center(struct multimin *data);
void multimin_find_best(struct multimin *data);
void sync(struct multimin *data,int s);
void multimin_shrink(struct multimin *data);
double contract(struct multimin *data,double mul);
double expand(struct multimin *data,double mul);
double reflect(struct multimin *data,double mul);
int simplex_iterate(struct simulation *sim,struct multimin *data);
void simplex_copy_ans(struct multimin *data);
void newton_start(struct simulation *sim,struct multimin *data);

#endif
