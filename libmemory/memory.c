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

/** @file memory.c
@brief get/free memory
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include <circuit.h>
#include "memory.h"
#include "ray_fun.h"
#include "newton_tricks.h"
#include "shape.h"

void device_alloc_traps(struct device *in)
{
	struct dimensions *dim=&in->ns.dim;

	//printf("hello %d %d %d %d \n",dim->xlen,dim->ylen,dim->zlen,dim->srh_bands);

	malloc_srh_bands(dim, &(in->nt));
	//printf("1\n");
	malloc_srh_bands(dim, &(in->ntlast));

	malloc_srh_bands(dim, &(in->dnt));
	//printf("2\n");
	malloc_srh_bands(dim, &(in->srh_n_r1));
	malloc_srh_bands(dim, &(in->srh_n_r2));
	malloc_srh_bands(dim, &(in->srh_n_r3));
	malloc_srh_bands(dim, &(in->srh_n_r4));
	malloc_srh_bands(dim, &(in->dsrh_n_r1));
	//printf("3\n");
	malloc_srh_bands(dim, &(in->dsrh_n_r2));
	malloc_srh_bands(dim, &(in->dsrh_n_r3));
	malloc_srh_bands(dim, &(in->dsrh_n_r4));
	malloc_srh_bands(dim, &(in->Fnt));
	malloc_srh_bands(dim, &(in->ntb_save));
	//printf("4\n");
	malloc_srh_bands(dim, &(in->nt_r1));
	malloc_srh_bands(dim, &(in->nt_r2));
	malloc_srh_bands(dim, &(in->nt_r3));
	malloc_srh_bands(dim, &(in->nt_r4));

	malloc_srh_bands(dim, &(in->pt));
	malloc_srh_bands(dim, &(in->ptlast));
	//printf("5\n");
	//printf("hello %d %d %d %d \n",dim->xlen,dim->ylen,dim->zlen,dim->srh_bands);

	malloc_srh_bands(dim, &(in->dpt));
	malloc_srh_bands(dim, &(in->srh_p_r1));
	malloc_srh_bands(dim, &(in->srh_p_r2));
	malloc_srh_bands(dim, &(in->srh_p_r3));
	malloc_srh_bands(dim, &(in->srh_p_r4));
	malloc_srh_bands(dim, &(in->dsrh_p_r1));
	malloc_srh_bands(dim, &(in->dsrh_p_r2));
	//printf("6\n");
	malloc_srh_bands(dim, &(in->dsrh_p_r3));
	//printf("7\n");
	malloc_srh_bands(dim, &(in->dsrh_p_r4));
	//printf("8\n");
	malloc_srh_bands(dim, &(in->ptb_save));
	//printf("9\n");
	malloc_srh_bands(dim, &(in->Fpt));
	//printf("10\n");
	malloc_srh_bands(dim, &(in->pt_r1));
	malloc_srh_bands(dim, &(in->pt_r2));
	malloc_srh_bands(dim, &(in->pt_r3));
	malloc_srh_bands(dim, &(in->pt_r4));

	//getchar();

	newton_state_alloc_traps(&(in->ns),dim);

}


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
	free_zxy_gdouble(dim,&in->Tl);
	free_zxy_gdouble(dim,&in->Te);
	free_zxy_gdouble(dim,&in->Th);
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
	free_zxy_gdouble(dim,&in->kl);
	free_zxy_gdouble(dim,&in->ke);
	free_zxy_gdouble(dim,&in->kh);
	free_zxy_gdouble(dim,&in->Hl);
	free_zxy_gdouble(dim,&in->He);
	free_zxy_gdouble(dim,&in->Hh);
	free_zxy_gdouble(dim,&in->Habs);
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


}

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

	malloc_zxy_gdouble(dim,&(in->Habs));

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

	malloc_zxy_gdouble(dim,&(in->Tl));

	malloc_zxy_gdouble(dim,&(in->Te));

	malloc_zxy_gdouble(dim,&(in->Th));

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

	malloc_zxy_gdouble(dim,&(in->kl));

	malloc_zxy_gdouble(dim,&(in->ke));

	malloc_zxy_gdouble(dim,&(in->kh));

	malloc_zxy_gdouble(dim,&(in->Hl));

	malloc_zxy_gdouble(dim,&(in->He));

	malloc_zxy_gdouble(dim,&(in->Hh));

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

	//newton_state_alloc_mesh(&(in->ns),dim);

	malloc_zxy_gdouble(dim,&(ns->phi));
	malloc_zxy_gdouble(dim,&(ns->x));
	malloc_zxy_gdouble(dim,&(ns->xp));

	circuit_alloc_nodes_and_links(sim,&(in->cir));

}
