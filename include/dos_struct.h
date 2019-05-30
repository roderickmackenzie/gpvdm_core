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
gdouble edge;
gdouble m;
int dostype;
gdouble Nt;
gdouble Et;
gdouble Nt2;
gdouble Et2;
gdouble Eshift;
gdouble nstart;
gdouble nstop;
gdouble base1;
gdouble base2;
int npoints;
int expan_len;
gdouble expan_N[20];
gdouble expan_E[20];
gdouble mu;
gdouble ion_density;
gdouble ion_mobility;
gdouble tau0;
gdouble tau1;
gdouble Tstart;
gdouble Tstop;
gdouble Ngaus;
int Tsteps;
int traps;
gdouble dband;
gdouble detrap;
int srh_bands;
gdouble srh_start;

gdouble srh_sigman;
gdouble srh_sigmap;
gdouble srh_vth;
gdouble Nc;
gdouble Nv;
gdouble srh_cut;

gdouble del_start;
gdouble del_stop;
gdouble Eg;
gdouble epsilonr;
gdouble doping_start;
gdouble doping_stop;
gdouble Xi;
gdouble gaus_mull;

gdouble pl_fe_fh;
gdouble pl_trap;
gdouble pl_recom;
int pl_enabled;

int Esteps;
gdouble B;
};

struct dos
{
int used;
gdouble *x;
int xlen;
int tlen;
int srh_bands;
gdouble *t;
gdouble *srh_E;
gdouble *srh_den;
gdouble **c;
gdouble **w;
gdouble ***srh_r1;
gdouble ***srh_r2;
gdouble ***srh_r3;
gdouble ***srh_r4;
gdouble ***srh_c;
struct dosconfig config;
};


#endif
