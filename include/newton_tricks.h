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

/** @file newton_tricks.h
@brief functions to solver for an external voltage
*/

#ifndef newton_tricks_h
#define newton_tricks_h


void newton_externv_aux(struct simulation *sim,struct device *in,gdouble V,gdouble* i,gdouble* didv,gdouble* didphi,gdouble* didxil,gdouble* didxipl,gdouble* didphir,gdouble* didxir,gdouble* didxipr);
gdouble newton_externv(struct simulation *sim,struct device *in,gdouble Vtot,int usecap);
gdouble newton_externalv_simple(struct simulation *sim,struct device *in,gdouble V);
long double sim_externalv_ittr(struct simulation *sim,struct device *in,gdouble wantedv);

//perovskite functions
long double newton_externalv_simple_perovskite(struct simulation *sim,struct device *in,gdouble V);
long double newton_externv_perovskite(struct simulation *sim,struct device *in,gdouble Vtot,int usecap);
#endif
