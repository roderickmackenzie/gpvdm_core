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

/** @file ntricks.c
	@brief A collection of helper functions for the newton solver.
*/


#include <stdio.h>
#include <exp.h>
#include "sim.h"
#include "dump.h"
#include "ntricks.h"
#include "gui_hooks.h"
#include <plot.h>
#include <cal_path.h>
#include <thermal.h>
#include <contacts.h>
#include <dump.h>
#include <log.h>

static int unused __attribute__((unused));

void save_state(struct simulation *sim,struct device *in,gdouble to)
{
	printf_log(sim,"Dumping state\n");
	int z=0;
	int x=0;
	int y=0;

	int band;
	FILE *state;
	state=fopena(get_output_path(sim),"state.dat","w");

	fprintf(state,"%Le ",to);

	for (z=0;z<in->zmeshpoints;z++)
	{
		for (x=0;x<in->xmeshpoints;x++)
		{
			for (y=0;y<in->ymeshpoints;y++)
			{
				fprintf(state,"%Le %Le %Le ",in->phi[z][x][y],in->x[z][x][y],in->xp[z][x][y]);

				for (band=0;band<in->srh_bands;band++)
				{
					fprintf(state,"%Le %Le ",in->xt[z][x][y][band],in->xpt[z][x][y][band]);
				}

			}
		}
	}

	fclose(state);
}

int load_state(struct simulation *sim,struct device *in,gdouble voltage)
{
	printf_log(sim,"Load state\n");
	int z=0;
	int x=0;
	int y=0;

	int band;
	gdouble vtest;
	FILE *state;
	state=fopena(get_output_path(sim),"state.dat","r");
	if (!state)
	{
	printf_log(sim,"State not found\n");
	return FALSE;
	}

	unused=fscanf(state,"%Le",&vtest);
	printf_log(sim,"%Le %Le",voltage,vtest);
	if (vtest!=voltage)
	{
	printf_log(sim,"State not found\n");
	return FALSE;
	}
	printf_log(sim,"Loading state\n");

	contact_set_active_contact_voltage(sim,in,vtest);

	for (z=0;z<in->zmeshpoints;z++)
	{
		for (x=0;x<in->xmeshpoints;x++)
		{
			for (y=0;y<in->ymeshpoints;y++)
			{
				unused=fscanf(state,"%Le %Le %Le ",&(in->phi[z][x][y]),&(in->x[z][x][y]),&(in->xp[z][x][y]));

				for (band=0;band<in->srh_bands;band++)
				{
					unused=fscanf(state,"%Le %Le ",&(in->xt[z][x][y][band]),&(in->xpt[z][x][y][band]));
				}

			}
		}
	}
	fclose(state);
return TRUE;

}


