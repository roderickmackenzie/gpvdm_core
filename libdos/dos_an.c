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

/** @file dos_an.c
	@brief Deals with the analytical DoS.
*/

#include "dos_an.h"
#include "const.h"
#include "code_ctrl.h"
#include "util.h"
#include "inp.h"
#include <log.h>


static int unused __attribute__((unused));


void dos_an_load(struct simulation *sim,struct dos_an_data *in,char *name)
{
//<clean>
strcpy(in->file,name);
struct inp_file inp;
inp_init(sim,&inp);
inp_load(sim,&inp,name);
inp_reset_read(sim,&inp);
in->items=0;
char *text;
do
{
	text=inp_get_string(sim,&inp);
	if (strcmp(text,"#end")==0) break;
 	if (strcmp(text,"#ver")==0) break;

	text=inp_get_string(sim,&inp);
	in->type[in->items]=english_to_bin(sim,text);

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	in->enable[in->items]=english_to_bin(sim,text);

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	sscanf(text,"%le",&(in->a[in->items]));

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	sscanf(text,"%le",&(in->b[in->items]));

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	sscanf(text,"%le",&(in->c[in->items]));

	in->items++;
}while(text!=NULL);

inp_free(sim,&inp);
//</clean>
}

double dos_an_get_value(struct simulation *sim,struct dos_an_data *in,double E)
{
//<clean>
int i=0;
double rho=0.0;
for (i=0;i<in->items;i++)
{
	if (in->enable[i]==TRUE)
	{
		if (in->type[i]==0)
		{
			in->a[i]=fabs(in->a[i]);
			in->b[i]=fabs(in->b[i]);
			in->c[i]=fabs(in->c[i]);
			rho+=in->a[i]*exp(-pow(((in->b[i]+E)/(sqrt(2.0)*in->c[i])),2.0));
		}else
		if (in->type[i]==1)
		{
			double Etail=fabs(in->b[i]);
			double Ntrap=fabs(in->a[i]);

			if (Etail<2e-3)
			{
				Etail=2e-3;
			}else
			if (Etail>200e-3)
			{
				 Etail=200e-3;
			}

			if (Ntrap<1e7)
			{
				Ntrap=1e7;
			}

			rho+=Ntrap*exp(E/Etail);
		}
		else
		{
			printf_log(sim,"In file %s DoS shape %d not known %d %d\n",in->file,in->type[i],i,in->items);
			exit(0);
		}
	}

}

return rho;
//</clean>
}

