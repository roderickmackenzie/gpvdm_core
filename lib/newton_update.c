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

/** @file newton_update.c
	@brief Update the dependent variables after running the newton solver to calculate new solution parameters.
*/

#include "sim.h"
#include "dump.h"
#include <dos.h>

void update_y_array(struct simulation *sim,struct device *in,int z,int x)
{
int y=0;
int band=0;
struct newton_save_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

	for (y=0;y<dim->ymeshpoints;y++)
	{
		in->Fn[z][x][y]=ns->x[z][x][y]-ns->phi[z][x][y];
		in->Fp[z][x][y]= -ns->xp[z][x][y]-ns->phi[z][x][y];

		in->Ec[z][x][y]= -ns->phi[z][x][y]-in->Xi[z][x][y];
		in->Ev[z][x][y]= -ns->phi[z][x][y]-in->Xi[z][x][y]-in->Eg[z][x][y];

		in->dn[z][x][y]=get_dn_den(in,ns->x[z][x][y]+in->t[z][x][y],in->Te[z][x][y],in->imat[z][x][y]);
		in->n[z][x][y]=get_n_den(in,ns->x[z][x][y]+in->t[z][x][y],in->Te[z][x][y],in->imat[z][x][y]);
		in->dndphi[z][x][y]=get_dn_den(in,ns->x[z][x][y]+in->t[z][x][y],in->Te[z][x][y],in->imat[z][x][y]);
		in->dp[z][x][y]=get_dp_den(in,ns->xp[z][x][y]-in->tp[z][x][y],in->Th[z][x][y],in->imat[z][x][y]);
		in->p[z][x][y]=get_p_den(in,ns->xp[z][x][y]-in->tp[z][x][y],in->Th[z][x][y],in->imat[z][x][y]);
		in->dpdphi[z][x][y]= -get_dp_den(in,ns->xp[z][x][y]-in->tp[z][x][y],in->Th[z][x][y],in->imat[z][x][y]);

		in->wn[z][x][y]=get_n_w(in,ns->x[z][x][y]+in->t[z][x][y],in->Te[z][x][y],in->imat[z][x][y]);
		in->wp[z][x][y]=get_p_w(in,ns->xp[z][x][y]-in->tp[z][x][y],in->Th[z][x][y],in->imat[z][x][y]);

		//in->mun[z][x][y]=get_n_mu(in,in->imat[z][x][y]);
		//in->mup[z][x][y]=get_p_mu(in,in->imat[z][x][y]);

		if (in->ntrapnewton)
		{
			in->nt_all[z][x][y]=0.0;
			for (band=0;band<dim->srh_bands;band++)
			{
				in->Fnt[z][x][y][band]=ns->xt[z][x][y][band]-ns->phi[z][x][y];

				in->srh_n_r1[z][x][y][band]=get_n_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_1,in->imat[z][x][y]);
				in->srh_n_r2[z][x][y][band]=get_n_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_2,in->imat[z][x][y]);
				in->srh_n_r3[z][x][y][band]=get_n_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_3,in->imat[z][x][y]);
				in->srh_n_r4[z][x][y][band]=get_n_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_4,in->imat[z][x][y]);
				in->dsrh_n_r1[z][x][y][band]=get_dn_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_1,in->imat[z][x][y]);
				in->dsrh_n_r2[z][x][y][band]=get_dn_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_2,in->imat[z][x][y]);
				in->dsrh_n_r3[z][x][y][band]=get_dn_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_3,in->imat[z][x][y]);
				in->dsrh_n_r4[z][x][y][band]=get_dn_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,srh_4,in->imat[z][x][y]);

				in->nt[z][x][y][band]=get_n_pop_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,in->imat[z][x][y]);
				in->dnt[z][x][y][band]=get_dn_pop_srh(sim,in,ns->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,in->imat[z][x][y]);
				in->nt_all[z][x][y]+=in->nt[z][x][y][band];

			}
		}

		if (in->ptrapnewton)
		{
			in->pt_all[z][x][y]=0.0;
			for (band=0;band<dim->srh_bands;band++)
			{
				in->Fpt[z][x][y][band]= -ns->xpt[z][x][y][band]-ns->phi[z][x][y];

				in->srh_p_r1[z][x][y][band]=get_p_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_1,in->imat[z][x][y]);
				in->srh_p_r2[z][x][y][band]=get_p_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_2,in->imat[z][x][y]);
				in->srh_p_r3[z][x][y][band]=get_p_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_3,in->imat[z][x][y]);
				in->srh_p_r4[z][x][y][band]=get_p_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_4,in->imat[z][x][y]);
				in->dsrh_p_r1[z][x][y][band]=get_dp_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_1,in->imat[z][x][y]);
				in->dsrh_p_r2[z][x][y][band]=get_dp_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_2,in->imat[z][x][y]);
				in->dsrh_p_r3[z][x][y][band]=get_dp_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_3,in->imat[z][x][y]);
				in->dsrh_p_r4[z][x][y][band]=get_dp_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,srh_4,in->imat[z][x][y]);

				in->pt[z][x][y][band]=get_p_pop_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,in->imat[z][x][y]);
				in->dpt[z][x][y][band]=get_dp_pop_srh(sim,in,ns->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,in->imat[z][x][y]);
				in->pt_all[z][x][y]+=in->pt[z][x][y][band];
			}
		}

		}

}


void init_mat_arrays(struct device *in)
{
int x=0;
int y=0;
int z=0;
struct newton_save_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

	for (z=0;z<dim->zmeshpoints;z++)
	{

		for (x=0;x<dim->xmeshpoints;x++)
		{

			for (y=0;y<dim->ymeshpoints;y++)
			{
				in->Tl[z][x][y]=in->Tll+dim->ymesh[y]*(in->Tlr-in->Tll)/in->ylen;
				in->Te[z][x][y]=in->Tll+dim->ymesh[y]*(in->Tlr-in->Tll)/in->ylen;
				in->Th[z][x][y]=in->Tll+dim->ymesh[y]*(in->Tlr-in->Tll)/in->ylen;
				in->ex[z][x][y]=0.0;
				in->Hex[z][x][y]=0.0;
				//if ((i>in->ymeshpoints/2)&&(i<in->ymeshpoints/2+10)) in->Hex[z][x][y]=1e9;
				in->epsilonr[z][x][y]=get_dos_epsilonr(in,in->imat[z][x][y]);

				in->Eg[z][x][y]=get_dos_Eg(in,in->imat[z][x][y]);

				in->muion[z][x][y]=get_dos_ion_mobility(in,in->imat[z][x][y]);
				in->Nion[z][x][y]=get_dos_ion_density(in,in->imat[z][x][y]);

				//printf("%d %d %d %Lf %d\n",z,x,y,in->Eg[z][x][y],in->imat[z][x][y]);
				//getchar();

				in->B[z][x][y]=get_dos_B(in,in->imat[z][x][y]);
				in->Dex[z][x][y]=0.0;//get_mat_param(&(in->mat.l[in->imat[z][x][y]]),mat_Dex);

				in->Xi[z][x][y]=get_dos_Xi(in,in->imat[z][x][y]);

				in->Ec[z][x][y]= -in->Xi[z][x][y];

				in->Ev[z][x][y]= -in->Xi[z][x][y]-in->Eg[z][x][y];


				in->Nc[z][x][y]=get_Nc_free(in,in->imat[z][x][y]);

				in->Nv[z][x][y]=get_Nv_free(in,in->imat[z][x][y]);

				in->mun[z][x][y]=get_n_mu(in,in->imat[z][x][y]);
				in->mup[z][x][y]=get_p_mu(in,in->imat[z][x][y]);

				in->kf[z][x][y]=0.0;//get_mat_param(&(in->mat.l[in->imat[z][x][y]]),mat_kf);
				in->kl[z][x][y]=in->thermal_kl;//get_mat_param(&(in->mat.l[in->imat[z][x][y]]),mat_kl);
				in->ke[z][x][y]=get_n_mu(in,in->imat[z][x][y]);
				in->kh[z][x][y]=get_p_mu(in,in->imat[z][x][y]);

				in->Hl[z][x][y]=0.0;
				in->He[z][x][y]=0.0;
				in->Hh[z][x][y]=0.0;
				in->Habs[z][x][y]=0.0;

				in->t[z][x][y]=in->Xi[z][x][y];
				in->tp[z][x][y]=in->Xi[z][x][y]+in->Eg[z][x][y];

				in->tt[z][x][y]=in->Xi[z][x][y];
				in->tpt[z][x][y]=in->Xi[z][x][y]+in->Eg[z][x][y];

			}
		}
	}
}

