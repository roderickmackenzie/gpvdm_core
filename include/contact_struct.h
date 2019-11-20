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

/** @file contact_struct.h
	@brief Definition of the contacts.
*/


#ifndef contact_struct_h
#define contact_struct_h

#include <shape_struct.h>

struct contact
{
	char name[100];
	int position;
	int active;
	gdouble voltage;
	long double voltage_want;
	gdouble voltage_last;
	gdouble store;
	long double np;
	int charge_type;
	long double area;

	char shape_file_name[100];
	struct shape shape;
	long double Rcontact;
	int type;
	long double ve0;
	long double vh0;
};

#endif
