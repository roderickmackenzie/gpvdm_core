//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/** @file dos_an.c
	@brief Deals with the analytical DoS.
*/

#include "dos_an.h"
#include "gpvdm_const.h"
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
	in->a[in->items]=fabs(in->a[in->items]);

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	sscanf(text,"%le",&(in->b[in->items]));
	in->b[in->items]=fabs(in->b[in->items]);

	text=inp_get_string(sim,&inp);
	text=inp_get_string(sim,&inp);
	sscanf(text,"%le",&(in->c[in->items]));
	in->c[in->items]=fabs(in->c[in->items]);

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

