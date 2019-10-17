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

/** @file dos_struct.h
	@brief Hold information about the DoS.
*/


#ifndef dos_struct_h
#define dos_struct_h

struct dosconfig
{
	char dos_name[20];
	char analytical_dos_file_name[20];
	int dos_number;
	long double m;
	int dostype;
	long double Nt;
	long double Et;
	long double nstart;
	long double nstop;
	int npoints;
	long double mu;
	long double ion_density;
	long double ion_mobility;
	long double Tstart;
	long double Tstop;
	int Tsteps;
	int traps;
	long double dband;
	long double detrap;
	int srh_bands;
	long double srh_start;

	long double srh_sigman;
	long double srh_sigmap; 
	long double srh_vth;
	long double Nc;
	long double Nv;
	long double srh_cut;

	long double del_start;
	long double del_stop;
	long double Eg;
	long double Xi;
	long double epsilonr;
	long double doping_start;
	long double doping_stop;

	int Esteps;
	long double B;
};

struct dos
{
	int used;
	long double *x;
	int xlen;
	int tlen;
	int srh_bands;
	long double *t;
	long double *srh_E;
	long double *srh_den;
	long double **c;
	long double **w;
	long double ***srh_r1;
	long double ***srh_r2;
	long double ***srh_r3;
	long double ***srh_r4;
	long double ***srh_c;
	struct dosconfig config;
};


#endif
