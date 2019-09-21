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

/** @file hash.c
	@brief Hashing function for fast lookup in tables, but the arrays are now linear so you don't need to hash.
*/

#include <stdio.h>
#include "sim.h"
#include "inp.h"
#include "cal_path.h"
#include "newton_tricks.h"

void hash_dir(struct simulation *sim)
{
	printf("here\n");
	struct list files;
	inp_listdir(sim, get_input_path(sim),&files);
	list_dump(&files);
	inp_list_free(&files);
}

