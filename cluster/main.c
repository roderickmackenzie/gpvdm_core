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

/** @file main.c
@brief main function
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "state.h"

struct state sim;

struct state *get_sim()
{
	return &sim;
}

int main(int argc, char *argv[])
{
packet_init_mutex();
printf("Clustering code\n");
state_init(&sim);
encrypt_load(&sim);

if (strcmp(argv[1],"--head")==0)
{
	sim.state=HEAD;
	head(&sim);
}

if (strcmp(argv[1],"--node")==0)
{
	sim.state=NODE;
	node(&sim);
}

}
