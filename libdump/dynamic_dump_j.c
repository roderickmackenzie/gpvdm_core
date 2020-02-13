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

/** @file dynamic_dump_j.c
@brief Dumps dynamic current density
*/

#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <i.h>
#include <exp.h>
#include <dos.h>
#include "sim.h"
#include "dump.h"
#include "dat_file.h"
#include "dynamic_store.h"
#include "memory.h"
#include "contacts.h"
#include <lang.h>
#include <cal_path.h>


static int unused __attribute__((unused));

void dump_dynamic_save_j(struct simulation *sim,struct device *in,char *outputpath,struct dynamic_store *store)
{
	int i;
	int sub=TRUE;
	char temp[200];
	struct dat_file buf;
	struct istruct one;

	buffer_init(&buf);



	if (get_dump_status(sim,dump_norm_time_to_one)==TRUE)
	{
		buf.norm_x_axis=TRUE;
	}

	if (get_dump_status(sim,dump_norm_y_axis)==TRUE)
	{
		buf.norm_y_axis=TRUE;
	}

	char out_dir[1000];
	join_path(2, out_dir,outputpath,"dynamic");
	struct stat st = {0};

	buffer_add_dir(sim,out_dir);

	char outpath[200];

	//////////////////////////////////////////J
	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->J_y0_n));
	sprintf(buf.title,"%s",_("Electron current density left contact (J_y0_n)"));
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->J_y0_n).x, (store->J_y0_n).data, (store->J_y0_n).len);
	buffer_dump_path(sim,out_dir,"J_y0_n.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->J_y0_p));
	sprintf(buf.title,"%s",_("Hole current density left contact (J_y0_p)"));
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->J_y0_p).x, (store->J_y0_p).data, (store->J_y0_p).len);
	buffer_dump_path(sim,out_dir,"J_y0_p.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->J_y1_n));
	sprintf(buf.title,"%s",_("Electron current density right contact (J_y1_n)"));
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->J_y1_n).x, (store->J_y1_n).data, (store->J_y1_n).len);
	buffer_dump_path(sim,out_dir,"J_y1_n.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->J_y1_p));
	sprintf(buf.title,"%s",_("Hole current density right contact (J_y1_p)"));
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->J_y1_p).x, (store->J_y1_p).data, (store->J_y1_p).len);
	buffer_dump_path(sim,out_dir,"J_y1_p.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jn_drift));
	sprintf(buf.title,"%s",_("Electron drift current"));
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jn_drift).x, (store->dynamic_jn_drift).data, (store->dynamic_jn_drift).len);
	buffer_dump_path(sim,out_dir,"jn_drift.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jn_diffusion));
	strcpy(buf.title,_("Electron diffusion current"));
	strcpy(buf.type,"xy");
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jn_diffusion).x, (store->dynamic_jn_diffusion).data, (store->dynamic_jn_diffusion).len);
	buffer_dump_path(sim,out_dir,"jn_diffusion.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jp_drift));
	sprintf(buf.title,"%s",_("Hole drift current"));
	strcpy(buf.type,"xy");
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jp_drift).x, (store->dynamic_jp_drift).data, (store->dynamic_jp_drift).len);
	buffer_dump_path(sim,out_dir,"jp_drift.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jp_diffusion));
	sprintf(buf.title,"%s",_("Hole diffusion current"));
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jp_diffusion).x, (store->dynamic_jp_diffusion).data, (store->dynamic_jp_diffusion).len);
	buffer_dump_path(sim,out_dir,"jp_diffusion.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jn));
	sprintf(buf.title,"%s",_("Jn at contacts"));
	strcpy(buf.type,"xy");
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jn).x, (store->dynamic_jn).data, (store->dynamic_jn).len);
	buffer_dump_path(sim,out_dir,"jn_contacts.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jp));
	strcpy(buf.title,_("Jp at contacts"));
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.y_units,"\\mu s");
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->dynamic_jp).x, (store->dynamic_jp).data, (store->dynamic_jp).len);
	buffer_dump_path(sim,out_dir,"jp_contacts.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,(store->jn_avg).x, (store->jn_avg).data, (store->jn_avg).len);
	buffer_dump_path(sim,out_dir,"jn_avg.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,(store->jp_avg).x, (store->jp_avg).data, (store->jp_avg).len);
	buffer_dump_path(sim,out_dir,"jp_avg.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->iout));
	strcpy(buf.title,_("External Current"));
	strcpy(buf.data_label,_("Current"));
	strcpy(buf.data_units,"Amps");
	buffer_add_info(sim,&buf);
	buffer_add_xy_data(sim,&buf,(store->iout).x, (store->iout).data, (store->iout).len);
	join_path(3, outpath,outputpath,"dynamic","i.dat");
	buffer_free(&buf);

	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,(store->jnout_mid).x, (store->jnout_mid).data, (store->jnout_mid).len);
	buffer_dump_path(sim,out_dir,"jn_mid.dat",&buf);
	buffer_free(&buf);


	inter_init(sim,&one);
	inter_copy(&one,&(store->jnout_mid),TRUE);
	inter_deriv(&one,&(store->jnout_mid));

	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,one.x, one.data, one.len);
	buffer_dump_path(sim,out_dir,"djn.dat",&buf);
	buffer_free(&buf);

	inter_free(&one);

	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,(store->jpout_mid).x, (store->jpout_mid).data, (store->jpout_mid).len);
	buffer_dump_path(sim,out_dir,"jp_mid.dat",&buf);
	buffer_free(&buf);

	inter_copy(&one,&(store->jpout_mid),TRUE);
	inter_deriv(&one,&(store->jpout_mid));
	buffer_malloc(&buf);
	buffer_add_xy_data(sim,&buf,one.x, one.data, one.len);
	buffer_dump_path(sim,out_dir,"djp.dat",&buf);
	buffer_free(&buf);
	inter_free(&one);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jp_drift));
	sprintf(buf.title,"%s + %s",_("Hole drift current"),_(" Hole diffusion current"));
	strcpy(buf.data_label,"Hole current density");
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);

	for (i=0;i<(store->dynamic_jp_drift).len;i++)
	{
		sprintf(temp,"%Le %Le\n",(store->dynamic_jp_drift).x[i],(store->dynamic_jp_drift).data[i]+(store->dynamic_jp_diffusion).data[i]);
		buffer_add_string(&buf,temp);
	}
	buffer_dump_path(sim,out_dir,"jp_drift_plus_diffusion.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->dynamic_jn_drift));
	sprintf(buf.title,"%s + %s",_("Electron drift current"),_("Electron diffusion current"));
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	for (i=0;i<(store->dynamic_jn_drift).len;i++)
	{
		sprintf(temp,"%Le %Le\n",(store->dynamic_jn_drift).x[i],(store->dynamic_jn_drift).data[i]+(store->dynamic_jn_diffusion).data[i]);
		buffer_add_string(&buf,temp);
	}
	buffer_dump_path(sim,out_dir,"jn_drift_plus_diffusion.dat",&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	buf.data_mul=1.0;
	dynamic_info_to_buf(sim,&buf, in,&(store->jout));
	strcpy(buf.title,_("Current density at contacts"));
	strcpy(buf.data_label,_("Current density"));
	strcpy(buf.data_units,"A m^{-2}");
	buffer_add_info(sim,&buf);
	if (sub==TRUE)
	{
		inter_sub_long_double(&(store->jout),(store->jout).data[0]);
		inter_mul(&(store->jout),-1.0);
	}
	buffer_add_xy_data(sim,&buf,(store->jout).x, (store->jout).data, (store->jout).len);
	buffer_dump_path(sim,out_dir,"j.dat",&buf);
	buffer_free(&buf);

}

