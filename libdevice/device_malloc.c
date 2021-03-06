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


/** @file device_malloc.c
	@brief Malloc for the device structure.
*/

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <device.h>
#include <string.h>
#include <dump.h>
#include <mesh.h>
#include <ray_fun.h>
#include <newton_tricks.h>
#include <memory.h>
#include <circuit.h>
#include <shape.h>
#include <lang.h>
#include <util.h>
#include <heat_fun.h>
#include <device_fun.h>

static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));

void device_get_memory(struct simulation *sim,struct device *in)
{
	in->odes = 0;
	struct newton_state *ns=(&in->ns);
	struct dimensions *dim=&in->ns.dim;

	in->obj=malloc(sizeof(struct object)*1000);

	if (dim->ylen==0)
	{
		return;
	}

	if ((dim->ylen<1)||(dim->xlen<1)||(dim->zlen<1))
	{
		ewe(sim,"%s\n",_("I can't allocate a device with less than 1 mesh point."));
	}

	if ((dim->ylen>50000)||(dim->xlen>50000)||(dim->zlen>50000))
	{
		ewe(sim,"%s\n",_("You are asking me to simulate a device with more than 50000 mesh points, although I could do this I am not going to because it seems a bad idea to me."));
	}



	//1d

	//2d
	malloc_zx_gdouble(dim,&(in->Vapplied_y0));
	malloc_zx_gdouble(dim,&(in->Vapplied_y1));

	malloc_zy_long_double(dim,&(in->Vapplied_x0));
	malloc_zy_long_double(dim,&(in->Vapplied_x1));

	malloc_zx_gdouble(dim,&(in->Jn_y0));
	malloc_zx_gdouble(dim,&(in->Jn_y1));
	malloc_zx_gdouble(dim,&(in->Jp_y0));
	malloc_zx_gdouble(dim,&(in->Jp_y1));

	malloc_zx_int(dim,&(in->n_contact_y0));
	malloc_zx_int(dim,&(in->n_contact_y1));
	malloc_zy_int(dim,&(in->n_contact_x0));
	malloc_zy_int(dim,&(in->n_contact_x1));

	malloc_zx_int(dim,&(in->passivate_y0));
	malloc_zx_int(dim,&(in->passivate_y1));
	malloc_zy_int(dim,&(in->passivate_x0));
	malloc_zy_int(dim,&(in->passivate_x1));

	malloc_zx_gdouble(dim,&(in->electrons_y0));
	malloc_zx_gdouble(dim,&(in->holes_y0));
	malloc_zx_gdouble(dim,&(in->electrons_y1));
	malloc_zx_gdouble(dim,&(in->holes_y1));

	malloc_zy_long_double(dim,&(in->electrons_x0));
	malloc_zy_long_double(dim,&(in->holes_x0));
	malloc_zy_long_double(dim,&(in->electrons_x1));
	malloc_zy_long_double(dim,&(in->holes_x1));


	malloc_zx_gdouble(dim,&(in->Fi0_y0));
	malloc_zx_gdouble(dim,&(in->Fi0_y1));
	malloc_zy_long_double(dim,&(in->Fi0_x0));
	malloc_zy_long_double(dim,&(in->Fi0_x1));

	malloc_zx_gdouble(dim,&(in->V_y0));
	malloc_zx_gdouble(dim,&(in->V_y1));
	malloc_zy_long_double(dim,&(in->V_x0));
	malloc_zy_long_double(dim,&(in->V_x1));


	//3d
	malloc_zxy_gdouble(dim,&(in->nf_save));

	malloc_zxy_gdouble(dim,&(in->pf_save));

	malloc_zxy_gdouble(dim,&(in->nt_save));

	malloc_zxy_gdouble(dim,&(in->pt_save));

	malloc_zxy_gdouble(dim,&(in->nfequlib));

	malloc_zxy_gdouble(dim,&(in->pfequlib));

	malloc_zxy_gdouble(dim,&(in->ntequlib));

	malloc_zxy_gdouble(dim,&(in->ptequlib));

	malloc_zxy_gdouble(dim,&(in->B));

	malloc_zxy_gdouble(dim,&(in->Nad));

	malloc_zxy_gdouble(dim,&(in->n));

	malloc_zxy_gdouble(dim,&(in->p));

	malloc_zxy_gdouble(dim,&(in->dn));

	malloc_zxy_gdouble(dim,&(in->dp));

	malloc_zxy_gdouble(dim,&(in->dndphi));

	malloc_zxy_gdouble(dim,&(in->dpdphi));

	malloc_zxy_gdouble(dim,&(in->Eg));

	malloc_zxy_gdouble(dim,&(in->Fn));

	malloc_zxy_gdouble(dim,&(in->Fp));

	malloc_zxy_gdouble(dim,&(in->Xi));

	malloc_zxy_gdouble(dim,&(in->Ev));

	malloc_zxy_gdouble(dim,&(in->Ec));

	malloc_zxy_gdouble(dim,&(in->mun));

	malloc_zxy_gdouble(dim,&(in->mup));

	malloc_zxy_gdouble(dim,&(in->muion));

	malloc_zxy_gdouble(dim,&(in->Nion));

	malloc_zxy_gdouble(dim,&(in->Nion_last));

	malloc_zxy_gdouble(dim,&(in->Dn));

	malloc_zxy_gdouble(dim,&(in->Dp));

	malloc_zxy_gdouble(dim,&(in->Nc));

	malloc_zxy_gdouble(dim,&(in->Nv));

	malloc_zxy_gdouble(dim,&(in->G));

	malloc_zxy_gdouble(dim,&(in->Gn));

	malloc_zxy_gdouble(dim,&(in->Photon_gen));

	malloc_zxy_gdouble(dim,&(in->Gp));


	malloc_zxy_gdouble(dim,&(in->Fi));

	malloc_zxy_gdouble(dim,&(in->Jn));

	malloc_zxy_gdouble(dim,&(in->Jp));

	malloc_zxy_gdouble(dim,&(in->Jn_x));

	malloc_zxy_gdouble(dim,&(in->Jp_x));

	malloc_zxy_gdouble(dim,&(in->Jn_drift));

	malloc_zxy_gdouble(dim,&(in->Jn_diffusion));

	malloc_zxy_gdouble(dim,&(in->Jp_drift));

	malloc_zxy_gdouble(dim,&(in->Jp_diffusion));


	malloc_zxy_gdouble(dim,&(in->t));

	malloc_zxy_gdouble(dim,&(in->tp));

	malloc_zxy_gdouble(dim,&(in->kf));

	malloc_zxy_gdouble(dim,&(in->kd));

	malloc_zxy_gdouble(dim,&(in->kr));

	malloc_zxy_gdouble(dim,&(in->Rfree));

	malloc_zxy_gdouble(dim,&(in->Rn));
	malloc_zxy_gdouble(dim,&(in->Rp));

	malloc_zxy_gdouble(dim,&(in->Rn_srh));
	malloc_zxy_gdouble(dim,&(in->Rp_srh));

	malloc_zxy_gdouble(dim,&(in->Rnet));

	malloc_zxy_gdouble(dim,&(in->ex));

	malloc_zxy_gdouble(dim,&(in->Dex));

	malloc_zxy_gdouble(dim,&(in->Hex));

	malloc_zxy_gdouble(dim,&(in->epsilonr));
	malloc_zxy_gdouble(dim,&(in->epsilonr_e0));



	malloc_zxy_gdouble(dim,&(in->nlast));

	malloc_zxy_gdouble(dim,&(in->plast));

	malloc_zxy_gdouble(dim,&(in->wn));

	malloc_zxy_gdouble(dim,&(in->wp));

	malloc_zxy_gdouble(dim,&(in->nt_all));

	malloc_zxy_gdouble(dim,&(in->phi_save));

	malloc_zxy_gdouble(dim,&(in->tt));


	malloc_zxy_gdouble(dim,&(in->pt_all));

	malloc_zxy_gdouble(dim,&(in->tpt));

	malloc_zxy_gdouble(dim,&(in->Rbi_k));

	malloc_zxy_gdouble(dim,&(in->nrelax));

	malloc_zxy_gdouble(dim,&(in->ntrap_to_p));

	malloc_zxy_gdouble(dim,&(in->prelax));

	malloc_zxy_gdouble(dim,&(in->ptrap_to_n));

	malloc_zxy_gdouble(dim,&(in->n_orig));

	malloc_zxy_gdouble(dim,&(in->p_orig));

	malloc_zxy_gdouble(dim,&(in->n_orig_f));

	malloc_zxy_gdouble(dim,&(in->p_orig_f));

	malloc_zxy_gdouble(dim,&(in->n_orig_t));

	malloc_zxy_gdouble(dim,&(in->p_orig_t));

	malloc_3d_int(dim,&(in->imat));
	malloc_3d_int(dim,&(in->imat_epitaxy));
	malloc_3d_int(dim,&(in->mask));

	dim_alloc_zx_epitaxy(&(in->dim_epitaxy),in);
	malloc_zx_epitaxy_int(&(in->dim_epitaxy),&(in->mask_epitaxy));

	newton_state_alloc_mesh(ns,dim, FALSE);
	//newton_state_alloc(ns,dim);
	//newton_state_alloc_mesh(&(in->ns),dim);

	//malloc_zxy_gdouble(dim,&(ns->phi));
	//malloc_zxy_gdouble(dim,&(ns->x));
	//malloc_zxy_gdouble(dim,&(ns->xp));

	circuit_alloc_nodes_and_links(sim,&(in->cir));

	//Thermal
	malloc_zxy_gdouble(dim,&(in->Tl));
	malloc_zxy_gdouble(dim,&(in->Te));
	malloc_zxy_gdouble(dim,&(in->Th));

	malloc_zxy_gdouble(dim,&(in->ke));
	malloc_zxy_gdouble(dim,&(in->kh));

	malloc_zxy_gdouble(dim,&(in->Hl));
	malloc_zxy_gdouble(dim,&(in->H_recombination));
	malloc_zxy_gdouble(dim,&(in->H_joule));

	malloc_zxy_gdouble(dim,&(in->He));
	malloc_zxy_gdouble(dim,&(in->Hh));

}

