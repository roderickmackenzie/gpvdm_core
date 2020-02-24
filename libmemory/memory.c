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

/** @file memory.c
@brief get/free memory
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include <circuit.h>
#include "memory.h"
#include "ray_fun.h"
#include "newton_tricks.h"
#include "shape.h"

void device_alloc_traps(struct device *in)
{
	struct dimensions *dim=&in->ns.dim;

	//printf("hello %d %d %d %d \n",dim->xlen,dim->ylen,dim->zlen,dim->srh_bands);

	malloc_srh_bands(dim, &(in->nt));
	//printf("1\n");
	malloc_srh_bands(dim, &(in->ntlast));

	malloc_srh_bands(dim, &(in->dnt));
	//printf("2\n");
	malloc_srh_bands(dim, &(in->srh_n_r1));
	malloc_srh_bands(dim, &(in->srh_n_r2));
	malloc_srh_bands(dim, &(in->srh_n_r3));
	malloc_srh_bands(dim, &(in->srh_n_r4));
	malloc_srh_bands(dim, &(in->dsrh_n_r1));
	//printf("3\n");
	malloc_srh_bands(dim, &(in->dsrh_n_r2));
	malloc_srh_bands(dim, &(in->dsrh_n_r3));
	malloc_srh_bands(dim, &(in->dsrh_n_r4));
	malloc_srh_bands(dim, &(in->Fnt));
	malloc_srh_bands(dim, &(in->ntb_save));
	//printf("4\n");
	malloc_srh_bands(dim, &(in->nt_r1));
	malloc_srh_bands(dim, &(in->nt_r2));
	malloc_srh_bands(dim, &(in->nt_r3));
	malloc_srh_bands(dim, &(in->nt_r4));

	malloc_srh_bands(dim, &(in->pt));
	malloc_srh_bands(dim, &(in->ptlast));
	//printf("5\n");
	//printf("hello %d %d %d %d \n",dim->xlen,dim->ylen,dim->zlen,dim->srh_bands);

	malloc_srh_bands(dim, &(in->dpt));
	malloc_srh_bands(dim, &(in->srh_p_r1));
	malloc_srh_bands(dim, &(in->srh_p_r2));
	malloc_srh_bands(dim, &(in->srh_p_r3));
	malloc_srh_bands(dim, &(in->srh_p_r4));
	malloc_srh_bands(dim, &(in->dsrh_p_r1));
	malloc_srh_bands(dim, &(in->dsrh_p_r2));
	//printf("6\n");
	malloc_srh_bands(dim, &(in->dsrh_p_r3));
	//printf("7\n");
	malloc_srh_bands(dim, &(in->dsrh_p_r4));
	//printf("8\n");
	malloc_srh_bands(dim, &(in->ptb_save));
	//printf("9\n");
	malloc_srh_bands(dim, &(in->Fpt));
	//printf("10\n");
	malloc_srh_bands(dim, &(in->pt_r1));
	malloc_srh_bands(dim, &(in->pt_r2));
	malloc_srh_bands(dim, &(in->pt_r3));
	malloc_srh_bands(dim, &(in->pt_r4));

	//getchar();

	newton_state_alloc_traps(&(in->ns),dim);

}


