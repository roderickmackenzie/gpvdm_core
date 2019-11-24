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
#include <complex_solver.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include "memory.h"
#include "ray_fun.h"
#include "newton_tricks.h"

void device_alloc_traps(struct device *in)
{
	struct dimensions *dim=&in->ns.dim;

	//printf("hello %d %d %d %d \n",dim->xmeshpoints,dim->ymeshpoints,dim->zmeshpoints,dim->srh_bands);

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
	//printf("hello %d %d %d %d \n",dim->xmeshpoints,dim->ymeshpoints,dim->zmeshpoints,dim->srh_bands);

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

	newton_save_state_alloc_traps(&(in->ns),dim);

}


void device_free(struct simulation *sim,struct device *in)
{
	struct dimensions *dim=&in->ns.dim;

	if (dim->ymeshpoints==0)
	{
		return;
	}

	//2d
	free_zx_gdouble(dim,&in->Vapplied_r);
	free_zx_gdouble(dim,&in->Vapplied_l);
	free_zx_gdouble(dim,&in->Jnleft);
	free_zx_gdouble(dim,&in->Jnright);
	free_zx_gdouble(dim,&in->Jpleft);
	free_zx_gdouble(dim,&in->Jpright);
	free_zx_int(dim,&in->n_contact_r);
	free_zx_int(dim,&in->n_contact_l);
	free_zx_int(dim,&in->passivate_r);
	free_zx_int(dim,&in->passivate_l);
	free_zx_gdouble(dim,&in->l_electrons);
	free_zx_gdouble(dim,&in->l_holes);
	free_zx_gdouble(dim,&in->r_electrons);
	free_zx_gdouble(dim,&in->r_holes);
	free_zx_gdouble(dim,&in->Fi0_top);
	free_zx_gdouble(dim,&in->Vl);
	free_zx_gdouble(dim,&in->Vr);

	//3d
	free_3d_gdouble(dim,&in->B);
	free_3d_gdouble(dim,&in->Nad);
	free_3d_gdouble(dim,&in->n);
	free_3d_gdouble(dim,&in->p);
	free_3d_gdouble(dim,&in->dn);
	free_3d_gdouble(dim,&in->dp);
	free_3d_gdouble(dim,&in->dndphi);
	free_3d_gdouble(dim,&in->dpdphi);
	free_3d_gdouble(dim,&in->Eg);
	free_3d_gdouble(dim,&in->Xi);
	free_3d_gdouble(dim,&in->Ev);
	free_3d_gdouble(dim,&in->Ec);
	free_3d_gdouble(dim,&in->mun);
	free_3d_gdouble(dim,&in->mup);
	free_3d_gdouble(dim,&in->muion);
	free_3d_gdouble(dim,&in->Nion);
	free_3d_gdouble(dim,&in->Nion_last);
	free_3d_gdouble(dim,&in->Dn);
	free_3d_gdouble(dim,&in->Dp);
	free_3d_gdouble(dim,&in->Fn);
	free_3d_gdouble(dim,&in->Fp);

	free_3d_gdouble(dim,&in->Nc);
	free_3d_gdouble(dim,&in->Nv);
	free_3d_gdouble(dim,&in->G);
	free_3d_gdouble(dim,&in->Gn);
	free_3d_gdouble(dim,&in->Gp);
	free_3d_gdouble(dim,&in->Photon_gen);
	free_3d_gdouble(dim,&in->Tl);
	free_3d_gdouble(dim,&in->Te);
	free_3d_gdouble(dim,&in->Th);
	free_3d_gdouble(dim,&in->Fi);
	free_3d_gdouble(dim,&in->Jn);
	free_3d_gdouble(dim,&in->Jp);
	free_3d_gdouble(dim,&in->Jn_x);
	free_3d_gdouble(dim,&in->Jp_x);
	free_3d_gdouble(dim,&in->Jn_drift);
	free_3d_gdouble(dim,&in->Jn_diffusion);
	free_3d_gdouble(dim,&in->Jp_drift);
	free_3d_gdouble(dim,&in->Jp_diffusion);
	free_3d_gdouble(dim,&in->t);
	free_3d_gdouble(dim,&in->tp);
	free_3d_gdouble(dim,&in->ex);
	free_3d_gdouble(dim,&in->Dex);
	free_3d_gdouble(dim,&in->Hex);
	free_3d_gdouble(dim,&in->epsilonr);

	free_3d_gdouble(dim,&in->kf);
	free_3d_gdouble(dim,&in->kd);
	free_3d_gdouble(dim,&in->kr);
	free_3d_gdouble(dim,&in->Rfree);
	free_3d_gdouble(dim,&in->Rn);
	free_3d_gdouble(dim,&in->Rp);
	free_3d_gdouble(dim,&in->Rn_srh);
	free_3d_gdouble(dim,&in->Rp_srh);
	free_3d_gdouble(dim,&in->Rnet);
	free_3d_gdouble(dim,&in->kl);
	free_3d_gdouble(dim,&in->ke);
	free_3d_gdouble(dim,&in->kh);
	free_3d_gdouble(dim,&in->Hl);
	free_3d_gdouble(dim,&in->He);
	free_3d_gdouble(dim,&in->Hh);
	free_3d_gdouble(dim,&in->Habs);
	free_3d_gdouble(dim,&in->nlast);
	free_3d_gdouble(dim,&in->plast);

	free_3d_gdouble(dim,&in->wn);
	free_3d_gdouble(dim,&in->wp);

	free_3d_gdouble(dim,&in->nt_all);

	free_3d_gdouble(dim,&in->tt);
	free_3d_gdouble(dim,&in->Rbi_k);

	free_3d_gdouble(dim,&in->pt_all);


	free_3d_gdouble(dim,&in->tpt);

	free_3d_gdouble(dim,&in->nf_save);
	free_3d_gdouble(dim,&in->pf_save);
	free_3d_gdouble(dim,&in->nt_save);
	free_3d_gdouble(dim,&in->pt_save);

	free_3d_gdouble(dim,&in->nfequlib);
	free_3d_gdouble(dim,&in->pfequlib);
	free_3d_gdouble(dim,&in->ntequlib);
	free_3d_gdouble(dim,&in->ptequlib);

	free_3d_gdouble(dim,&in->nrelax);
	free_3d_gdouble(dim,&in->ntrap_to_p);
	free_3d_gdouble(dim,&in->prelax);
	free_3d_gdouble(dim,&in->ptrap_to_n);

	free_3d_gdouble(dim,&in->n_orig);
	free_3d_gdouble(dim,&in->p_orig);
	free_3d_gdouble(dim,&in->n_orig_f);
	free_3d_gdouble(dim,&in->p_orig_f);
	free_3d_gdouble(dim,&in->n_orig_t);
	free_3d_gdouble(dim,&in->p_orig_t);

	free_3d_gdouble(dim,&in->phi_save);

	free_3d_int(dim,in->imat);
	free_3d_int(dim,in->imat_epitaxy);

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

	newton_save_state_free(&(in->ns));
	//dim_free(&(in->dim_max));
	//Free epitaxy

	//Free solvers
	solver_free(sim);
	complex_solver_free(sim);
	printf_log(sim,"%s %i %s\n", _("Solved"), in->odes, _("Equations"));

}

void device_get_memory(struct simulation *sim,struct device *in)
{
	in->odes = 0;
	in->Ti = NULL;
	in->Tj = NULL;
	in->Tx = NULL;
	in->b = NULL;
	in->Tdebug = NULL;


	struct dimensions *dim=&in->ns.dim;

	//dim_alloc(dim);
	if (dim->ymeshpoints==0)
	{
		return;
	}

	if ((dim->ymeshpoints<1)||(dim->xmeshpoints<1)||(dim->zmeshpoints<1))
	{
		ewe(sim,"%s\n",_("I can't allocate a device with less than 1 mesh point."));
	}

	if ((dim->ymeshpoints>50000)||(dim->xmeshpoints>50000)||(dim->zmeshpoints>50000))
	{
		ewe(sim,"%s\n",_("You are asking me to simulate a device with more than 50000 mesh points, although I could do this I am not going to because it seems a bad idea to me."));
	}



	//1d


	//2d
	malloc_zx_gdouble(dim,&(in->Vapplied_r));
	malloc_zx_gdouble(dim,&(in->Vapplied_l));
	malloc_zx_gdouble(dim,&(in->Jnleft));
	malloc_zx_gdouble(dim,&(in->Jnright));
	malloc_zx_gdouble(dim,&(in->Jpleft));
	malloc_zx_gdouble(dim,&(in->Jpright));
	malloc_zx_int(dim,&(in->n_contact_r));
	malloc_zx_int(dim,&(in->n_contact_l));
	malloc_zx_int(dim,&(in->passivate_r));
	malloc_zx_int(dim,&(in->passivate_l));

	malloc_zx_gdouble(dim,&(in->l_electrons));
	malloc_zx_gdouble(dim,&(in->l_holes));
	malloc_zx_gdouble(dim,&(in->r_electrons));
	malloc_zx_gdouble(dim,&(in->r_holes));
	malloc_zx_gdouble(dim,&(in->Fi0_top));
	malloc_zx_gdouble(dim,&(in->Vl));
	malloc_zx_gdouble(dim,&(in->Vr));


	//3d
	zxy_malloc_gdouble(dim,&(in->nf_save));

	zxy_malloc_gdouble(dim,&(in->pf_save));

	zxy_malloc_gdouble(dim,&(in->nt_save));

	zxy_malloc_gdouble(dim,&(in->pt_save));

	zxy_malloc_gdouble(dim,&(in->nfequlib));

	zxy_malloc_gdouble(dim,&(in->pfequlib));

	zxy_malloc_gdouble(dim,&(in->ntequlib));

	zxy_malloc_gdouble(dim,&(in->ptequlib));

	zxy_malloc_gdouble(dim,&(in->Habs));

	zxy_malloc_gdouble(dim,&(in->B));

	zxy_malloc_gdouble(dim,&(in->Nad));

	zxy_malloc_gdouble(dim,&(in->n));

	zxy_malloc_gdouble(dim,&(in->p));

	zxy_malloc_gdouble(dim,&(in->dn));

	zxy_malloc_gdouble(dim,&(in->dp));

	zxy_malloc_gdouble(dim,&(in->dndphi));

	zxy_malloc_gdouble(dim,&(in->dpdphi));

	zxy_malloc_gdouble(dim,&(in->Eg));

	zxy_malloc_gdouble(dim,&(in->Fn));

	zxy_malloc_gdouble(dim,&(in->Fp));

	zxy_malloc_gdouble(dim,&(in->Xi));

	zxy_malloc_gdouble(dim,&(in->Ev));

	zxy_malloc_gdouble(dim,&(in->Ec));

	zxy_malloc_gdouble(dim,&(in->mun));

	zxy_malloc_gdouble(dim,&(in->mup));

	zxy_malloc_gdouble(dim,&(in->muion));

	zxy_malloc_gdouble(dim,&(in->Nion));

	zxy_malloc_gdouble(dim,&(in->Nion_last));

	zxy_malloc_gdouble(dim,&(in->Dn));

	zxy_malloc_gdouble(dim,&(in->Dp));

	zxy_malloc_gdouble(dim,&(in->Nc));

	zxy_malloc_gdouble(dim,&(in->Nv));

	zxy_malloc_gdouble(dim,&(in->G));

	zxy_malloc_gdouble(dim,&(in->Gn));

	zxy_malloc_gdouble(dim,&(in->Photon_gen));

	zxy_malloc_gdouble(dim,&(in->Gp));

	zxy_malloc_gdouble(dim,&(in->Tl));

	zxy_malloc_gdouble(dim,&(in->Te));

	zxy_malloc_gdouble(dim,&(in->Th));

	zxy_malloc_gdouble(dim,&(in->Fi));

	zxy_malloc_gdouble(dim,&(in->Jn));

	zxy_malloc_gdouble(dim,&(in->Jp));

	zxy_malloc_gdouble(dim,&(in->Jn_x));

	zxy_malloc_gdouble(dim,&(in->Jp_x));

	zxy_malloc_gdouble(dim,&(in->Jn_drift));

	zxy_malloc_gdouble(dim,&(in->Jn_diffusion));

	zxy_malloc_gdouble(dim,&(in->Jp_drift));

	zxy_malloc_gdouble(dim,&(in->Jp_diffusion));


	zxy_malloc_gdouble(dim,&(in->t));

	zxy_malloc_gdouble(dim,&(in->tp));

	zxy_malloc_gdouble(dim,&(in->kf));

	zxy_malloc_gdouble(dim,&(in->kd));

	zxy_malloc_gdouble(dim,&(in->kr));

	zxy_malloc_gdouble(dim,&(in->Rfree));

	zxy_malloc_gdouble(dim,&(in->Rn));
	zxy_malloc_gdouble(dim,&(in->Rp));

	zxy_malloc_gdouble(dim,&(in->Rn_srh));
	zxy_malloc_gdouble(dim,&(in->Rp_srh));

	zxy_malloc_gdouble(dim,&(in->Rnet));

	zxy_malloc_gdouble(dim,&(in->ex));

	zxy_malloc_gdouble(dim,&(in->Dex));

	zxy_malloc_gdouble(dim,&(in->Hex));

	zxy_malloc_gdouble(dim,&(in->epsilonr));

	zxy_malloc_gdouble(dim,&(in->kl));

	zxy_malloc_gdouble(dim,&(in->ke));

	zxy_malloc_gdouble(dim,&(in->kh));

	zxy_malloc_gdouble(dim,&(in->Hl));

	zxy_malloc_gdouble(dim,&(in->He));

	zxy_malloc_gdouble(dim,&(in->Hh));

	zxy_malloc_gdouble(dim,&(in->nlast));

	zxy_malloc_gdouble(dim,&(in->plast));

	zxy_malloc_gdouble(dim,&(in->wn));

	zxy_malloc_gdouble(dim,&(in->wp));

	zxy_malloc_gdouble(dim,&(in->nt_all));

	zxy_malloc_gdouble(dim,&(in->phi_save));

	zxy_malloc_gdouble(dim,&(in->tt));


	zxy_malloc_gdouble(dim,&(in->pt_all));

	zxy_malloc_gdouble(dim,&(in->tpt));

	zxy_malloc_gdouble(dim,&(in->Rbi_k));

	zxy_malloc_gdouble(dim,&(in->nrelax));

	zxy_malloc_gdouble(dim,&(in->ntrap_to_p));

	zxy_malloc_gdouble(dim,&(in->prelax));

	zxy_malloc_gdouble(dim,&(in->ptrap_to_n));

	zxy_malloc_gdouble(dim,&(in->n_orig));

	zxy_malloc_gdouble(dim,&(in->p_orig));

	zxy_malloc_gdouble(dim,&(in->n_orig_f));

	zxy_malloc_gdouble(dim,&(in->p_orig_f));

	zxy_malloc_gdouble(dim,&(in->n_orig_t));

	zxy_malloc_gdouble(dim,&(in->p_orig_t));

	malloc_3d_int(dim,&(in->imat));
	malloc_3d_int(dim,&(in->imat_epitaxy));

	newton_save_state_alloc_mesh(&(in->ns),dim);


}
