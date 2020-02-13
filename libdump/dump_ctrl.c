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

/** @file dump_ctrl.c
@brief set/get what is allowed to be written to disk, writing to disk is bad.
*/

#include <dump_ctrl.h>
#include <gpvdm_const.h>
#include <sim.h>

int get_dump_status(struct simulation *sim,int a)
{
	return sim->dump_array[a];
}

void set_dump_status(struct simulation *sim,int name, int a)
{
if ((sim->dump_array[dump_lock]==FALSE)||(name==dump_lock)) sim->dump_array[name]=a;
}

