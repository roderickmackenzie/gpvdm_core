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
	struct dimensions *dim=&in->dim;

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
	struct dimensions *dim=&in->dim;

	if (dim->ymeshpoints==0)
	{
		return;
	}

	//2d
	free_zx_gdouble(dim,in->Vapplied_r);
	free_zx_gdouble(dim,in->Vapplied_l);
	free_zx_gdouble(dim,in->Jnleft);
	free_zx_gdouble(dim,in->Jnright);
	free_zx_gdouble(dim,in->Jpleft);
	free_zx_gdouble(dim,in->Jpright);
	free_zx_int(dim,in->n_contact_r);
	free_zx_int(dim,in->n_contact_l);
	free_zx_int(dim,in->passivate_r);
	free_zx_int(dim,in->passivate_l);
	free_zx_gdouble(dim,in->l_electrons);
	free_zx_gdouble(dim,in->l_holes);
	free_zx_gdouble(dim,in->r_electrons);
	free_zx_gdouble(dim,in->r_holes);
	free_zx_gdouble(dim,in->Fi0_top);
	free_zx_gdouble(dim,in->Vl);
	free_zx_gdouble(dim,in->Vr);

	//3d
	free_3d_gdouble(dim,in->B);
	free_3d_gdouble(dim,in->Nad);
	free_3d_gdouble(dim,in->n);
	free_3d_gdouble(dim,in->p);
	free_3d_gdouble(dim,in->dn);
	free_3d_gdouble(dim,in->dp);
	free_3d_gdouble(dim,in->dndphi);
	free_3d_gdouble(dim,in->dpdphi);
	free_3d_gdouble(dim,in->Eg);
	free_3d_gdouble(dim,in->Xi);
	free_3d_gdouble(dim,in->Ev);
	free_3d_gdouble(dim,in->Ec);
	free_3d_gdouble(dim,in->mun);
	free_3d_gdouble(dim,in->mup);
	free_3d_gdouble(dim,in->muion);
	free_3d_gdouble(dim,in->Nion);
	free_3d_gdouble(dim,in->Nion_last);
	free_3d_gdouble(dim,in->Dn);
	free_3d_gdouble(dim,in->Dp);
	free_3d_gdouble(dim,in->Fn);
	free_3d_gdouble(dim,in->Fp);

	free_3d_gdouble(dim,in->Nc);
	free_3d_gdouble(dim,in->Nv);
	free_3d_gdouble(dim,in->G);
	free_3d_gdouble(dim,in->Gn);
	free_3d_gdouble(dim,in->Gp);
	free_3d_gdouble(dim,in->Photon_gen);
	free_3d_gdouble(dim,in->Tl);
	free_3d_gdouble(dim,in->Te);
	free_3d_gdouble(dim,in->Th);
	free_3d_gdouble(dim,in->R);
	free_3d_gdouble(dim,in->Fi);
	free_3d_gdouble(dim,in->Jn);
	free_3d_gdouble(dim,in->Jp);
	free_3d_gdouble(dim,in->Jn_x);
	free_3d_gdouble(dim,in->Jp_x);
	free_3d_gdouble(dim,in->Jn_drift);
	free_3d_gdouble(dim,in->Jn_diffusion);
	free_3d_gdouble(dim,in->Jp_drift);
	free_3d_gdouble(dim,in->Jp_diffusion);
	free_3d_gdouble(dim,in->t);
	free_3d_gdouble(dim,in->tp);
	free_3d_gdouble(dim,in->ex);
	free_3d_gdouble(dim,in->Dex);
	free_3d_gdouble(dim,in->Hex);
	free_3d_gdouble(dim,in->epsilonr);

	free_3d_gdouble(dim,in->kf);
	free_3d_gdouble(dim,in->kd);
	free_3d_gdouble(dim,in->kr);
	free_3d_gdouble(dim,in->Rfree);
	free_3d_gdouble(dim,in->Rn);
	free_3d_gdouble(dim,in->Rp);
	free_3d_gdouble(dim,in->Rn_srh);
	free_3d_gdouble(dim,in->Rp_srh);
	free_3d_gdouble(dim,in->kl);
	free_3d_gdouble(dim,in->ke);
	free_3d_gdouble(dim,in->kh);
	free_3d_gdouble(dim,in->Hl);
	free_3d_gdouble(dim,in->He);
	free_3d_gdouble(dim,in->Hh);
	free_3d_gdouble(dim,in->Habs);
	free_3d_gdouble(dim,in->nlast);
	free_3d_gdouble(dim,in->plast);

	free_3d_gdouble(dim,in->wn);
	free_3d_gdouble(dim,in->wp);

	free_3d_gdouble(dim,in->nt_all);

	free_3d_gdouble(dim,in->tt);
	free_3d_gdouble(dim,in->Rbi_k);

	free_3d_gdouble(dim,in->pt_all);


	free_3d_gdouble(dim,in->tpt);

	free_3d_gdouble(dim,in->nf_save);
	free_3d_gdouble(dim,in->pf_save);
	free_3d_gdouble(dim,in->nt_save);
	free_3d_gdouble(dim,in->pt_save);

	free_3d_gdouble(dim,in->nfequlib);
	free_3d_gdouble(dim,in->pfequlib);
	free_3d_gdouble(dim,in->ntequlib);
	free_3d_gdouble(dim,in->ptequlib);

	free_3d_gdouble(dim,in->nrelax);
	free_3d_gdouble(dim,in->ntrap_to_p);
	free_3d_gdouble(dim,in->prelax);
	free_3d_gdouble(dim,in->ptrap_to_n);

	free_3d_gdouble(dim,in->n_orig);
	free_3d_gdouble(dim,in->p_orig);
	free_3d_gdouble(dim,in->n_orig_f);
	free_3d_gdouble(dim,in->p_orig_f);
	free_3d_gdouble(dim,in->n_orig_t);
	free_3d_gdouble(dim,in->p_orig_t);

	free_3d_gdouble(dim,in->phi_save);

	free_3d_int(dim,in->imat);
	free_3d_int(dim,in->imat_epitaxy);

	//traps
	free_srh_bands(dim, in->nt);
	free_srh_bands(dim, in->dnt);
	free_srh_bands(dim, in->srh_n_r1);
	free_srh_bands(dim, in->srh_n_r2);
	free_srh_bands(dim, in->srh_n_r3);
	free_srh_bands(dim, in->srh_n_r4);
	free_srh_bands(dim, in->dsrh_n_r1);
	free_srh_bands(dim, in->dsrh_n_r2);
	free_srh_bands(dim, in->dsrh_n_r3);
	free_srh_bands(dim, in->dsrh_n_r4);
	free_srh_bands(dim, in->Fnt);
	free_srh_bands(dim, in->ntb_save);

	free_srh_bands(dim, in->nt_r1);
	free_srh_bands(dim, in->nt_r2);
	free_srh_bands(dim, in->nt_r3);
	free_srh_bands(dim, in->nt_r4);

	free_srh_bands(dim, in->ntlast);

	free_srh_bands(dim, in->pt);
	free_srh_bands(dim, in->dpt);
	free_srh_bands(dim, in->srh_p_r1);
	free_srh_bands(dim, in->srh_p_r2);
	free_srh_bands(dim, in->srh_p_r3);
	free_srh_bands(dim, in->srh_p_r4);
	free_srh_bands(dim, in->dsrh_p_r1);
	free_srh_bands(dim, in->dsrh_p_r2);
	free_srh_bands(dim, in->dsrh_p_r3);
	free_srh_bands(dim, in->dsrh_p_r4);
	free_srh_bands(dim, in->Fpt);
	free_srh_bands(dim, in->ptb_save);

	free_srh_bands(dim, in->pt_r1);
	free_srh_bands(dim, in->pt_r2);
	free_srh_bands(dim, in->pt_r3);
	free_srh_bands(dim, in->pt_r4);

	free_srh_bands(dim, in->ptlast);

	newton_save_state_free(&(in->ns),&(in->dim));
	dim_free(&(in->dim));
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


	struct dimensions *dim=&in->dim;

	dim->xmeshpoints=in->mesh_data.meshdata_x.tot_points;
	dim->ymeshpoints=in->mesh_data.meshdata_y.tot_points;
	dim->zmeshpoints=in->mesh_data.meshdata_z.tot_points;

	dim_alloc(dim);
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
	malloc_3d_gdouble(dim,&(in->nf_save));

	malloc_3d_gdouble(dim,&(in->pf_save));

	malloc_3d_gdouble(dim,&(in->nt_save));

	malloc_3d_gdouble(dim,&(in->pt_save));

	malloc_3d_gdouble(dim,&(in->nfequlib));

	malloc_3d_gdouble(dim,&(in->pfequlib));

	malloc_3d_gdouble(dim,&(in->ntequlib));

	malloc_3d_gdouble(dim,&(in->ptequlib));

	malloc_3d_gdouble(dim,&(in->Habs));

	malloc_3d_gdouble(dim,&(in->B));

	malloc_3d_gdouble(dim,&(in->Nad));

	malloc_3d_gdouble(dim,&(in->n));

	malloc_3d_gdouble(dim,&(in->p));

	malloc_3d_gdouble(dim,&(in->dn));

	malloc_3d_gdouble(dim,&(in->dp));

	malloc_3d_gdouble(dim,&(in->dndphi));

	malloc_3d_gdouble(dim,&(in->dpdphi));

	malloc_3d_gdouble(dim,&(in->Eg));

	malloc_3d_gdouble(dim,&(in->Fn));

	malloc_3d_gdouble(dim,&(in->Fp));

	malloc_3d_gdouble(dim,&(in->Xi));

	malloc_3d_gdouble(dim,&(in->Ev));

	malloc_3d_gdouble(dim,&(in->Ec));

	malloc_3d_gdouble(dim,&(in->mun));

	malloc_3d_gdouble(dim,&(in->mup));

	malloc_3d_gdouble(dim,&(in->muion));

	malloc_3d_gdouble(dim,&(in->Nion));

	malloc_3d_gdouble(dim,&(in->Nion_last));

	malloc_3d_gdouble(dim,&(in->Dn));

	malloc_3d_gdouble(dim,&(in->Dp));

	malloc_3d_gdouble(dim,&(in->Nc));

	malloc_3d_gdouble(dim,&(in->Nv));

	malloc_3d_gdouble(dim,&(in->G));

	malloc_3d_gdouble(dim,&(in->Gn));

	malloc_3d_gdouble(dim,&(in->Photon_gen));

	malloc_3d_gdouble(dim,&(in->Gp));

	malloc_3d_gdouble(dim,&(in->Tl));

	malloc_3d_gdouble(dim,&(in->Te));

	malloc_3d_gdouble(dim,&(in->Th));

	malloc_3d_gdouble(dim,&(in->R));

	malloc_3d_gdouble(dim,&(in->Fi));

	malloc_3d_gdouble(dim,&(in->Jn));

	malloc_3d_gdouble(dim,&(in->Jp));

	malloc_3d_gdouble(dim,&(in->Jn_x));

	malloc_3d_gdouble(dim,&(in->Jp_x));

	malloc_3d_gdouble(dim,&(in->Jn_drift));

	malloc_3d_gdouble(dim,&(in->Jn_diffusion));

	malloc_3d_gdouble(dim,&(in->Jp_drift));

	malloc_3d_gdouble(dim,&(in->Jp_diffusion));


	malloc_3d_gdouble(dim,&(in->t));

	malloc_3d_gdouble(dim,&(in->tp));

	malloc_3d_gdouble(dim,&(in->kf));

	malloc_3d_gdouble(dim,&(in->kd));

	malloc_3d_gdouble(dim,&(in->kr));

	malloc_3d_gdouble(dim,&(in->Rfree));

	malloc_3d_gdouble(dim,&(in->Rn));
	malloc_3d_gdouble(dim,&(in->Rp));

	malloc_3d_gdouble(dim,&(in->Rn_srh));
	malloc_3d_gdouble(dim,&(in->Rp_srh));

	malloc_3d_gdouble(dim,&(in->ex));

	malloc_3d_gdouble(dim,&(in->Dex));

	malloc_3d_gdouble(dim,&(in->Hex));

	malloc_3d_gdouble(dim,&(in->epsilonr));

	malloc_3d_gdouble(dim,&(in->kl));

	malloc_3d_gdouble(dim,&(in->ke));

	malloc_3d_gdouble(dim,&(in->kh));

	malloc_3d_gdouble(dim,&(in->Hl));

	malloc_3d_gdouble(dim,&(in->He));

	malloc_3d_gdouble(dim,&(in->Hh));

	malloc_3d_gdouble(dim,&(in->nlast));

	malloc_3d_gdouble(dim,&(in->plast));

	malloc_3d_gdouble(dim,&(in->wn));

	malloc_3d_gdouble(dim,&(in->wp));

	malloc_3d_gdouble(dim,&(in->nt_all));

	malloc_3d_gdouble(dim,&(in->phi_save));

	malloc_3d_gdouble(dim,&(in->tt));


	malloc_3d_gdouble(dim,&(in->pt_all));

	malloc_3d_gdouble(dim,&(in->tpt));

	malloc_3d_gdouble(dim,&(in->Rbi_k));

	malloc_3d_gdouble(dim,&(in->nrelax));

	malloc_3d_gdouble(dim,&(in->ntrap_to_p));

	malloc_3d_gdouble(dim,&(in->prelax));

	malloc_3d_gdouble(dim,&(in->ptrap_to_n));

	malloc_3d_gdouble(dim,&(in->n_orig));

	malloc_3d_gdouble(dim,&(in->p_orig));

	malloc_3d_gdouble(dim,&(in->n_orig_f));

	malloc_3d_gdouble(dim,&(in->p_orig_f));

	malloc_3d_gdouble(dim,&(in->n_orig_t));

	malloc_3d_gdouble(dim,&(in->p_orig_t));

	malloc_3d_int(dim,&(in->imat));
	malloc_3d_int(dim,&(in->imat_epitaxy));

	newton_save_state_alloc_mesh(&(in->ns),&(in->dim));


}
