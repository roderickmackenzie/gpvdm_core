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

/** @file dump_contacts.c
@brief dump JV curves from the contacts.
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

#include "memory.h"
#include "contacts.h"
#include <lang.h>
#include <cal_path.h>
#include "contacts_vti_store.h"

static int unused __attribute__((unused));

void dump_contacts_init(struct simulation *sim,struct device *in,struct contacts_vti_store *store)
{
	if (in->ncontacts>2)
	{
		int i=0;

		for (i=0;i<in->ncontacts;i++)
		{
			inter_init(sim,&(store->J[i]));
		}
	}
}

void dump_contacts_save(struct simulation *sim,struct device *in,struct contacts_vti_store *store)
{
	char string[200];
	if (in->ncontacts>2)
	{
		int i;
		int ii;
		int sub=TRUE;
		char temp[200];
		struct dat_file buf;
		buffer_init(&buf);
		for (i=0;i<in->ncontacts;i++)
		{
			buffer_malloc(&buf);
			buf.y_mul=1.0;
			buf.data_mul=1.0;
			sprintf(buf.title,"%s (%s)",_("Voltage - Current density"), in->contacts[i].name);
			strcpy(buf.type,"xy");
			strcpy(buf.y_label,_("Voltage"));
			strcpy(buf.y_units,"V");
			strcpy(buf.data_label,_("Current"));
			strcpy(buf.data_units,"A m^{-2}");
			buf.logscale_y=0;
			buf.logscale_data=0;
			buf.x=1;
			buf.y=store->J[i].len;
			buf.z=1;
			buffer_add_info(sim,&buf);
			buffer_add_xy_data(sim,&buf,store->J[i].x, store->J[i].data, store->J[i].len);
			sprintf(temp,"jv_contact%d.dat",i);
			buffer_dump_path(sim,sim->output_path,temp,&buf);
			buffer_free(&buf);

			buffer_malloc(&buf);
			buf.y_mul=1.0;
			buf.data_mul=1.0;
			sprintf(buf.title,"%s (%s)",_("Voltage - Current"), in->contacts[i].name);
			strcpy(buf.type,"xy");
			strcpy(buf.y_label,_("Voltage"));
			strcpy(buf.y_units,"V");
			strcpy(buf.data_label,_("Current"));
			strcpy(buf.data_units,"A m^{-2}");
			buf.logscale_y=0;
			buf.logscale_data=0;
			buf.x=1;
			buf.y=store->J[i].len;
			buf.z=1;
			buffer_add_info(sim,&buf);

			for (ii=0;ii<store->J[i].len;ii++)
			{
				sprintf(string,"%Le %Le\n",store->J[i].x[ii],store->J[i].data[ii]*in->contacts[i].area);
				buffer_add_string(&buf,string);
			}

			sprintf(temp,"iv_contact%d.dat",i);
			buffer_dump_path(sim,sim->output_path,temp,&buf);
			buffer_free(&buf);

		}
	}
}

void dump_contacts_add_data(struct simulation *sim,struct device *in,struct contacts_vti_store *store)
{
	if (in->ncontacts>2)
	{
		int i=0;
		gdouble x_value=0.0;
		long double J=0.0;
		long double V_loss=0.0;
		for (i=0;i<in->ncontacts;i++)
		{
			J=contacts_get_J(in,i);
			x_value=contact_get_active_contact_voltage(sim,in);
			V_loss=in->contacts[i].area*J*10.0;//in->Rcontact;
	//		inter_append(&(store->v),x_value,contact_get_voltage(sim,in,i));
			inter_append(&(store->J[i]),x_value,contacts_get_J(in,i));
		}

	}

}

void dump_contacts_free(struct simulation *sim,struct device *in,struct contacts_vti_store *store)
{
	if (in->ncontacts>2)
	{
		int i=0;

		for (i=0;i<in->ncontacts;i++)
		{
			//inter_free(&(store->v[i]));
			inter_free(&(store->J[i]));
		}
	}
}

