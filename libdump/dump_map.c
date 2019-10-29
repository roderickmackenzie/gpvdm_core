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

/** @file dump_map.c
@brief dump a map charge carrier density map of the device, not really used any more due to python back end needs rewritign anyway
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dos.h>
#include "sim.h"
#include "dump.h"
#include "dat_file.h"
#include <cal_path.h>
#include <exp.h>
#include <lang.h>
#include <i.h>
#include <data2d.h>

static int unused __attribute__((unused));


void dump_device_map(struct simulation *sim,char* out_dir,struct device *in)
{
if (in->srh_bands==0) return;

if (in->dump_1d_slice_xpos>=in->xmeshpoints) return;
if (in->dump_1d_slice_zpos>=in->zmeshpoints) return;

int x=0;
int y=0;
struct dat_file buf;
char temp[1000];
int i=0;
long double Vexternal=get_equiv_V(sim,in);
int Epoints=0;
struct data2d map;
long double delta_E=dos_get_band_energy_n(in,1,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][0])-dos_get_band_energy_n(in,0,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][0]);
long double start=in->map_start;
long double stop=in->map_stop;
Epoints=(int)abs(stop-start)/delta_E;
int band=0;
long double E=0;
int vpos=0;
struct newton_save_state *ns=&(in->ns);

data2d_init(&map,in->ymeshpoints,Epoints);
data2d_init_y_mesh(&map,start, stop);


buffer_init(&buf);


//Electrons
buffer_malloc(&buf);
buf.y_mul=1e9;
buf.x_mul=1.0;
strcpy(buf.title,"Charge carrier density - position");
strcpy(buf.type,"heat");
strcpy(buf.x_label,_("Position"));
strcpy(buf.y_label,_("Energy"));
strcpy(buf.x_units,"nm");
strcpy(buf.y_units,"eV");

strcpy(buf.data_label,_("Charge density"));
strcpy(buf.data_units,"m^{-3} eV^{-1}");

buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_data=TRUE;

buf.x=Epoints;
buf.y=in->ymeshpoints;
buf.z=1;

buf.Vexternal=Vexternal;
buf.time=in->time;

buffer_add_info(sim,&buf);

data2d_set_value(&map,0.0);
for (band=0;band<in->srh_bands;band++)
{

	for (y=0;y<in->ymeshpoints;y++)
	{
		E=in->Ec[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]+dos_get_band_energy_n(in,band,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]);
		vpos=search(map.y_mesh,map.y_len,E);
		map.data[y][vpos]=in->nt[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y][band];
	}

}

for (y=0;y<Epoints;y++)
{

	for (x=0;x<in->ymeshpoints;x++)
	{
		sprintf(temp,"%Le %Le %Le\n",map.y_mesh[y],ns->ymesh[x],map.data[x][y]);
		buffer_add_string(&buf,temp);
	}

	buffer_add_string(&buf,"\n");
}


buffer_dump_path(sim,out_dir,"nt_map.dat",&buf);
buffer_free(&buf);





//Holes
buffer_malloc(&buf);

buf.y_mul=1e9;
buf.x_mul=1.0;
strcpy(buf.title,"Charge carrier density - position");
strcpy(buf.type,"heat");
strcpy(buf.x_label,_("Position"));
strcpy(buf.y_label,_("Energy"));
strcpy(buf.data_label,_("Charge density"));
strcpy(buf.x_units,"nm");
strcpy(buf.y_units,"eV");

strcpy(buf.data_label,_("Charge density"));
strcpy(buf.data_units,"m^{-3} eV^{-1}");

buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_data=TRUE;

buf.x=Epoints;
buf.y=in->ymeshpoints;
buf.z=1;

buf.Vexternal=Vexternal;
buf.time=in->time;

buffer_add_info(sim,&buf);

data2d_set_value(&map,0.0);
for (band=0;band<in->srh_bands;band++)
{

	for (y=0;y<in->ymeshpoints;y++)
	{
		E=in->Ev[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]-dos_get_band_energy_p(in,band,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]);
		vpos=search(map.y_mesh,map.y_len,E);
		map.data[y][vpos]=in->pt[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y][band];
	}

}

for (y=0;y<Epoints;y++)
{

	for (x=0;x<in->ymeshpoints;x++)
	{
		sprintf(temp,"%Le %Le %Le\n",map.y_mesh[y],ns->ymesh[x],map.data[x][y]);
		buffer_add_string(&buf,temp);
	}

	buffer_add_string(&buf,"\n");
}



buffer_dump_path(sim,out_dir,"pt_map.dat",&buf);

buffer_free(&buf);


//Charte total
buffer_malloc(&buf);
buf.y_mul=1e9;
buf.x_mul=1.0;
strcpy(buf.title,"Total charge - position");
strcpy(buf.type,"heat");
strcpy(buf.x_label,_("Position"));
strcpy(buf.y_label,_("Energy"));
strcpy(buf.x_units,"nm");
strcpy(buf.y_units,"eV");

strcpy(buf.data_label,_("Charge density"));
strcpy(buf.data_units,"m^{-3} eV^{-1}");

buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_data=TRUE;

buf.x=Epoints;
buf.y=in->ymeshpoints;
buf.z=1;

buf.Vexternal=Vexternal;
buf.time=in->time;

buffer_add_info(sim,&buf);


data2d_set_value(&map,0.0);
for (band=0;band<in->srh_bands;band++)
{

	for (y=0;y<in->ymeshpoints;y++)
	{
		E=in->Ev[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]-dos_get_band_energy_p(in,band,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]);
		vpos=search(map.y_mesh,map.y_len,E);
		map.data[y][vpos]=in->nt[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y][band];

		E=in->Ec[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]+dos_get_band_energy_n(in,band,in->imat[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y]);
		vpos=search(map.y_mesh,map.y_len,E);
		map.data[y][vpos]=in->pt[in->dump_1d_slice_zpos][in->dump_1d_slice_xpos][y][band];
	}

}

for (y=0;y<Epoints;y++)
{

	for (x=0;x<in->ymeshpoints;x++)
	{
		sprintf(temp,"%Le %Le %Le\n",map.y_mesh[y],ns->ymesh[x],map.data[x][y]);
		buffer_add_string(&buf,temp);
	}

	buffer_add_string(&buf,"\n");
}


buffer_dump_path(sim,out_dir,"npt_map.dat",&buf);
buffer_free(&buf);


data2d_free(&map);
}



