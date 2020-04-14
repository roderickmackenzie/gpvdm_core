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

/** @file contacts.c
@brief backend to handle complex contacts
*/

#include <string.h>
#include "epitaxy.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include <cal_path.h>
#include "contacts.h"
#include <log.h>
#include <dump.h>
#include <shape.h>
#include <shape_struct.h>

void contacts_load(struct simulation *sim,struct device *in)
{
	int i;
	struct inp_file inp;
	in->ncontacts=0;

	inp_init(sim,&inp);
	if (inp_load(sim, &inp , "contacts.inp")!=0)
	{
		ewe(sim,"Can't open the file contacts\n");
	}

	inp_check(sim,&inp,1.3);
	inp_reset_read(sim,&inp);
	inp_get_string(sim,&inp);
	sscanf(inp_get_string(sim,&inp),"%d",&(in->ncontacts));

	if (in->ncontacts>10)
	{
		ewe(sim,"Too many contacts\n");
	}

	int active=FALSE;
	long double ingress=0.0;

	struct shape *s;

	for (i=0;i<in->ncontacts;i++)
	{

		//inp_get_string(sim,&inp);	//name
		//strcpy(in->contacts[i].name,inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//position
		in->contacts[i].position=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//applied voltage type
		strcpy(in->contacts[i].applied_voltage_type,inp_get_string(sim,&inp));
		if (strcmp(in->contacts[i].applied_voltage_type,"change")==0)
		{
			in->contacts[i].active=TRUE;
		}else
		{
			in->contacts[i].active=FALSE;
		}

		if (strcmp(in->contacts[i].applied_voltage_type,"ground")==0)
		{
			in->contacts[i].ground=TRUE;
		}else
		{
			in->contacts[i].ground=FALSE;
		}


		inp_get_string(sim,&inp);	//applied_voltage
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].voltage_want));
		in->contacts[i].voltage=0.0;
		in->contacts[i].voltage_last=in->contacts[i].voltage;
		//printf("%Le\n",in->contacts[i].voltage_want);
		//getchar();

		inp_get_string(sim,&inp);	//np
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].np));
		in->contacts[i].np=fabs(in->contacts[i].np);

		inp_get_string(sim,&inp);	//charge_type
		in->contacts[i].charge_type=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//shape_file_name
		strcpy(in->contacts[i].shape_file_name,inp_get_string(sim,&inp));

		s=shape_load_file(sim,&(in->my_epitaxy),&(in->contacts[i].shape),in->contacts[i].shape_file_name);
		s->nx=1;
		s->nz=1;
		s->dz=in->zlen;
		strcpy(in->contacts[i].name,s->name);

		inp_get_string(sim,&inp);	//contact_resistance_sq
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].contact_resistance_sq));
		in->contacts[i].contact_resistance_sq=fabs(in->contacts[i].contact_resistance_sq);

		inp_get_string(sim,&inp);	//shunt_resistance_sq
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].shunt_resistance_sq));
		in->contacts[i].shunt_resistance_sq=fabs(in->contacts[i].shunt_resistance_sq);

		inp_get_string(sim,&inp);	//Physical model
		in->contacts[i].type=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//ve0
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].ve0));
		in->contacts[i].ve0=fabs(in->contacts[i].ve0);

		inp_get_string(sim,&inp);	//ve0
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].vh0));
		in->contacts[i].vh0=fabs(in->contacts[i].vh0);

		if (in->contacts[i].position==LEFT)
		{
			s->dx=ingress;
		}
	}

	char * ver = inp_get_string(sim,&inp);
	if (strcmp(ver,"#ver")!=0)
	{
			ewe(sim,"No #ver tag found in file\n");
	}

	inp_free(sim,&inp);

	contacts_update(sim,in);
	contact_set_flip_current(sim,in);

	contacts_cal_area(sim,in);

}

