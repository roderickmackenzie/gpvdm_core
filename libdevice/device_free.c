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


/** @file device_free.c
	@brief Free memory for the device structure.
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
#include <solver_interface.h>
#include <log.h>

static int unused __attribute__((unused));
static char* unused_pchar __attribute__((unused));


void device_free(struct simulation *sim,struct device *in)
{
	int i;
	struct dimensions *dim=&in->ns.dim;

	//world made from triangles
	for (i=0;i<in->objects;i++)
	{
		if (in->obj[i].n!=NULL)
		{
			ewe(sim,"There is still data in the array\n");
		}

		object_free(&(in->obj[i]));
	}

	free(in->obj);

	in->obj=NULL;
	in->objects=0;

	heat_free_memory(sim,&(in->thermal));

	shape_free(sim,&(in->big_box));

	if (dim->ylen==0)
	{
		return;
	}

	//2d

	free_zx_gdouble(dim,&in->Vapplied_y0);
	free_zx_gdouble(dim,&in->Vapplied_y1);

	free_zy_long_double(dim,&in->Vapplied_x0);
	free_zy_long_double(dim,&in->Vapplied_x1);

	free_zx_gdouble(dim,&in->Jn_y0);
	free_zx_gdouble(dim,&in->Jn_y1);
	free_zx_gdouble(dim,&in->Jp_y0);
	free_zx_gdouble(dim,&in->Jp_y1);

	free_zx_int(dim,&in->n_contact_y0);
	free_zx_int(dim,&in->n_contact_y1);
	free_zy_int(dim,&in->n_contact_x0);
	free_zy_int(dim,&in->n_contact_x1);


	free_zx_int(dim,&in->passivate_y0);
	free_zx_int(dim,&in->passivate_y1);
	free_zy_int(dim,&in->passivate_x0);
	free_zy_int(dim,&in->passivate_x1);

	free_zx_gdouble(dim,&in->electrons_y0);
	free_zx_gdouble(dim,&in->holes_y0);
	free_zx_gdouble(dim,&in->electrons_y1);
	free_zx_gdouble(dim,&in->holes_y1);

	free_zy_long_double(dim,&in->electrons_x0);
	free_zy_long_double(dim,&in->holes_x0);
	free_zy_long_double(dim,&in->electrons_x1);
	free_zy_long_double(dim,&in->holes_x1);

	free_zx_gdouble(dim,&in->Fi0_y0);
	free_zx_gdouble(dim,&in->Fi0_y1);
	free_zy_long_double(dim,&in->Fi0_x0);
	free_zy_long_double(dim,&in->Fi0_y1);

	free_zx_gdouble(dim,&in->V_y0);
	free_zx_gdouble(dim,&in->V_y1);
	free_zy_long_double(dim,&in->V_x0);
	free_zy_long_double(dim,&in->V_x1);

	//3d
	free_zxy_gdouble(dim,&in->B);
	free_zxy_gdouble(dim,&in->Nad);
	free_zxy_gdouble(dim,&in->n);
	free_zxy_gdouble(dim,&in->p);
	free_zxy_gdouble(dim,&in->dn);
	free_zxy_gdouble(dim,&in->dp);
	free_zxy_gdouble(dim,&in->dndphi);
	free_zxy_gdouble(dim,&in->dpdphi);
	free_zxy_gdouble(dim,&in->Eg);
	free_zxy_gdouble(dim,&in->Xi);
	free_zxy_gdouble(dim,&in->Ev);
	free_zxy_gdouble(dim,&in->Ec);
	free_zxy_gdouble(dim,&in->mun);
	free_zxy_gdouble(dim,&in->mup);
	free_zxy_gdouble(dim,&in->muion);
	free_zxy_gdouble(dim,&in->Nion);
	free_zxy_gdouble(dim,&in->Nion_last);
	free_zxy_gdouble(dim,&in->Dn);
	free_zxy_gdouble(dim,&in->Dp);
	free_zxy_gdouble(dim,&in->Fn);
	free_zxy_gdouble(dim,&in->Fp);

	free_zxy_gdouble(dim,&in->Nc);
	free_zxy_gdouble(dim,&in->Nv);
	free_zxy_gdouble(dim,&in->G);
	free_zxy_gdouble(dim,&in->Gn);
	free_zxy_gdouble(dim,&in->Gp);
	free_zxy_gdouble(dim,&in->Photon_gen);

	free_zxy_gdouble(dim,&in->Fi);
	free_zxy_gdouble(dim,&in->Jn);
	free_zxy_gdouble(dim,&in->Jp);
	free_zxy_gdouble(dim,&in->Jn_x);
	free_zxy_gdouble(dim,&in->Jp_x);
	free_zxy_gdouble(dim,&in->Jn_drift);
	free_zxy_gdouble(dim,&in->Jn_diffusion);
	free_zxy_gdouble(dim,&in->Jp_drift);
	free_zxy_gdouble(dim,&in->Jp_diffusion);
	free_zxy_gdouble(dim,&in->t);
	free_zxy_gdouble(dim,&in->tp);
	free_zxy_gdouble(dim,&in->ex);
	free_zxy_gdouble(dim,&in->Dex);
	free_zxy_gdouble(dim,&in->Hex);
	free_zxy_gdouble(dim,&in->epsilonr);

	free_zxy_gdouble(dim,&in->kf);
	free_zxy_gdouble(dim,&in->kd);
	free_zxy_gdouble(dim,&in->kr);
	free_zxy_gdouble(dim,&in->Rfree);
	free_zxy_gdouble(dim,&in->Rn);
	free_zxy_gdouble(dim,&in->Rp);
	free_zxy_gdouble(dim,&in->Rn_srh);
	free_zxy_gdouble(dim,&in->Rp_srh);
	free_zxy_gdouble(dim,&in->Rnet);
	free_zxy_gdouble(dim,&in->nlast);
	free_zxy_gdouble(dim,&in->plast);

	free_zxy_gdouble(dim,&in->wn);
	free_zxy_gdouble(dim,&in->wp);

	free_zxy_gdouble(dim,&in->nt_all);

	free_zxy_gdouble(dim,&in->tt);
	free_zxy_gdouble(dim,&in->Rbi_k);

	free_zxy_gdouble(dim,&in->pt_all);


	free_zxy_gdouble(dim,&in->tpt);

	free_zxy_gdouble(dim,&in->nf_save);
	free_zxy_gdouble(dim,&in->pf_save);
	free_zxy_gdouble(dim,&in->nt_save);
	free_zxy_gdouble(dim,&in->pt_save);

	free_zxy_gdouble(dim,&in->nfequlib);
	free_zxy_gdouble(dim,&in->pfequlib);
	free_zxy_gdouble(dim,&in->ntequlib);
	free_zxy_gdouble(dim,&in->ptequlib);

	free_zxy_gdouble(dim,&in->nrelax);
	free_zxy_gdouble(dim,&in->ntrap_to_p);
	free_zxy_gdouble(dim,&in->prelax);
	free_zxy_gdouble(dim,&in->ptrap_to_n);

	free_zxy_gdouble(dim,&in->n_orig);
	free_zxy_gdouble(dim,&in->p_orig);
	free_zxy_gdouble(dim,&in->n_orig_f);
	free_zxy_gdouble(dim,&in->p_orig_f);
	free_zxy_gdouble(dim,&in->n_orig_t);
	free_zxy_gdouble(dim,&in->p_orig_t);

	free_zxy_gdouble(dim,&in->phi_save);

	free_3d_int(dim,in->imat);
	free_3d_int(dim,in->imat_epitaxy);
	free_3d_int(dim,in->mask);

	//traps
	free_srh_bands(dim, &in->nt);
	free_srh_bands(dim, &in->dnt);
	free_srh_bands(dim, &in->srh_n_r1);
	free_srh_bands(dim, &in->srh_n_r2);
	free_srh_bands(dim, &in->srh_n_r3);
	free_srh_bands(dim, &in->srh_n_r4);
	free_srh_bands(dim, &in->dsrh_n_r1);
	free_srh_bands(dim, &in->dsrh_n_r2);
	free_srh_bands(dim, &in->dsrh_n_r3);
	free_srh_bands(dim, &in->dsrh_n_r4);
	free_srh_bands(dim, &in->Fnt);
	free_srh_bands(dim, &in->ntb_save);

	free_srh_bands(dim, &in->nt_r1);
	free_srh_bands(dim, &in->nt_r2);
	free_srh_bands(dim, &in->nt_r3);
	free_srh_bands(dim, &in->nt_r4);

	free_srh_bands(dim, &in->ntlast);

	free_srh_bands(dim, &in->pt);
	free_srh_bands(dim, &in->dpt);
	free_srh_bands(dim, &in->srh_p_r1);
	free_srh_bands(dim, &in->srh_p_r2);
	free_srh_bands(dim, &in->srh_p_r3);
	free_srh_bands(dim, &in->srh_p_r4);
	free_srh_bands(dim, &in->dsrh_p_r1);
	free_srh_bands(dim, &in->dsrh_p_r2);
	free_srh_bands(dim, &in->dsrh_p_r3);
	free_srh_bands(dim, &in->dsrh_p_r4);
	free_srh_bands(dim, &in->Fpt);
	free_srh_bands(dim, &in->ptb_save);

	free_srh_bands(dim, &in->pt_r1);
	free_srh_bands(dim, &in->pt_r2);
	free_srh_bands(dim, &in->pt_r3);
	free_srh_bands(dim, &in->pt_r4);

	free_srh_bands(dim, &in->ptlast);

	newton_state_free(&(in->ns));
	//dim_free(&(in->dim_max));
	//Free epitaxy

	//Free solvers
	solver_free(sim);
	printf_log(sim,"%s %i %s\n", _("Solved"), in->odes, _("Equations"));


	circuit_free(sim,&(in->cir));

	//Thermal
	free_zxy_gdouble(dim,&in->Tl);
	free_zxy_gdouble(dim,&in->Te);
	free_zxy_gdouble(dim,&in->Th);

	free_zxy_gdouble(dim,&in->Hl);
	free_zxy_gdouble(dim,&in->He);
	free_zxy_gdouble(dim,&in->Hh);

	free_zxy_gdouble(dim,&in->ke);
	free_zxy_gdouble(dim,&in->kh);

}


