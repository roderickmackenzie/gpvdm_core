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

/** @file newton.c
	@brief Standard newton solver.
*/


#include <string.h>
#include <log.h>
#include <gpvdm_const.h>
#include "newton.h"
#include <dll_export.h>
#include <util.h>
#include <exp.h>
#include <advmath.h>
#include <dump.h>
#include <cal_path.h>
#include <dos.h>
#include <sim.h>
#include <solver_interface.h>
#include <contacts.h>
#include <memory.h>
#include <newton_tricks.h>

static gdouble Jnl=0.0;
static gdouble Jnr=0.0;
static gdouble Jpl=0.0;
static gdouble Jpr=0.0;

static gdouble Dnl=0.0;
static gdouble Dnc=0.0;
static gdouble Dnr=0.0;
static gdouble Dpl=0.0;
static gdouble Dpc=0.0;
static gdouble Dpr=0.0;

static gdouble nl=0.0;
static gdouble nc=0.0;
static gdouble nr=0.0;

static gdouble pl=0.0;
static gdouble pc=0.0;
static gdouble pr=0.0;

static gdouble xil=0.0;
static gdouble xir=0.0;
static gdouble xipl=0.0;
static gdouble xipr=0.0;

static gdouble dJpdxipl=0.0;
static gdouble dJpdxipc=0.0;
static gdouble dJpdxipr=0.0;

static gdouble dnl=0.0;
static gdouble dnc=0.0;
static gdouble dnr=0.0;

static gdouble dpl=0.0;
static gdouble dpc=0.0;
static gdouble dpr=0.0;

static gdouble munl=0.0;
static gdouble munc=0.0;
static gdouble munr=0.0;

static gdouble mupl=0.0;
static gdouble mupc=0.0;
static gdouble mupr=0.0;


static gdouble wnl=0.0;
static gdouble wnc=0.0;
static gdouble wnr=0.0;

static gdouble wpl=0.0;
static gdouble wpc=0.0;
static gdouble wpr=0.0;

static gdouble dJdxil=0.0;
static gdouble dJdxic=0.0;
static gdouble dJdxir=0.0;

static gdouble dJdphil=0.0;
static gdouble dJdphic=0.0;
static gdouble dJdphir=0.0;

static gdouble dJpdphil=0.0;
static gdouble dJpdphic=0.0;
static gdouble dJpdphir=0.0;


static gdouble dphidxic=0.0;
static gdouble dphidxipc=0.0;


void update_solver_vars(struct simulation *sim,struct device *in,int z,int x,int clamp)
{
int i;
int band=0;
long double Vapplied=0.0;
long double clamp_temp=300.0;
struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;
struct matrix *mx=&in->mx;

gdouble update=0.0;

	for (i=0;i<dim->ylen;i++)
	{

		update=(gdouble)mx->b[i];
		if (clamp==TRUE)
		{
			ns->phi[z][x][i]+=update/(1.0+gfabs(update/in->electrical_clamp/(clamp_temp*kb/Q)));
		}else
		{
			ns->phi[z][x][i]+=update;
		}


		update=(gdouble)(mx->b[dim->ylen*(1)+i]);
		if (clamp==TRUE)
		{
			ns->x[z][x][i]+=update/(1.0+gfabs(update/in->electrical_clamp/(clamp_temp*kb/Q)));
		}else
		{
			ns->x[z][x][i]+=update;
		}


		update=(gdouble)(mx->b[dim->ylen*(1+1)+i]);
		if (clamp==TRUE)
		{
			ns->xp[z][x][i]+=update/(1.0+gfabs(update/in->electrical_clamp/(clamp_temp*kb/Q)));
		}else
		{
			ns->xp[z][x][i]+=update;

		}


		if (in->ntrapnewton==TRUE)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				update=(gdouble)(mx->b[dim->ylen*(1+1+1+band)+i]);
				if (clamp==TRUE)
				{
					ns->xt[z][x][i][band]+=update/(1.0+gfabs(update/in->electrical_clamp/(clamp_temp*kb/Q)));

				}else
				{
					ns->xt[z][x][i][band]+=update;
				}
			}
		}

		if (in->ptrapnewton==TRUE)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				update=(gdouble)(mx->b[dim->ylen*(1+1+1+dim->srh_bands+band)+i]);
				if (clamp==TRUE)
				{
					ns->xpt[z][x][i][band]+=update/(1.0+gfabs(update/in->electrical_clamp/(clamp_temp*kb/Q)));
				}else
				{
					ns->xpt[z][x][i][band]+=update;

				}
			}
		}


		}


update_y_array(sim,in,z,x);

}

void fill_matrix(struct simulation *sim,struct device *in,int z,int x)
{

struct matrix *mx=&(in->mx);
struct heat *thermal=&(in->thermal);

int band=0;
update_y_array(sim,in,z,x);

//FILE *file_j =fopen("myj.dat","w");
//getchar();
gdouble phil;
gdouble phic;
gdouble phir;
gdouble yl;
gdouble yc;
gdouble yr;
gdouble dyl;
gdouble dyr;
gdouble ddh=0.0;
//gdouble dh;
int pos=0;

gdouble Ecl=0.0;
gdouble Ecr=0.0;
gdouble Ecc=0.0;
gdouble Evl=0.0;
gdouble Evr=0.0;
gdouble Evc=0.0;

gdouble Tel=0.0;
//gdouble Tec=0.0;
gdouble Ter=0.0;

gdouble Thl=0.0;
//gdouble Thc=0.0;
gdouble Thr=0.0;

gdouble xnr;
gdouble tnr;
gdouble xnl;
gdouble tnl;

gdouble xpr;
gdouble tpr;
gdouble xpl;
gdouble tpl;


//gdouble exl;
//gdouble exr;
//gdouble exc;
//gdouble Dexl;
//gdouble Dexc;
//gdouble Dexr;
//gdouble R;

gdouble epr;
gdouble epc;
gdouble epl;

//gdouble G;
gdouble Gn;
gdouble Gp;
int i;
gdouble dJdxipc=0.0;
gdouble dJpdxic=0.0;

gdouble e0=0.0;
gdouble e1=0.0;

gdouble dphil=0.0;
gdouble dphic=0.0;
gdouble dphir=0.0;
gdouble deriv=0.0;

gdouble nlast=0.0;
gdouble plast=0.0;
gdouble dt=0.0;

//gdouble kll=0.0;
//gdouble klc=0.0;
//gdouble klr=0.0;

//gdouble kl0=0.0;
//gdouble kl1=0.0;

gdouble one=0.0;
//gdouble one0_l=0.0;
//gdouble one0_r=0.0;



gdouble Rtrapn=0.0;
gdouble Rtrapp=0.0;

gdouble dJndphil_leftl=0.0;
gdouble dJndphil_leftc=0.0;
gdouble dJpdphil_leftl=0.0;
gdouble dJpdphil_leftc=0.0;
gdouble dphil_left=0.0;
gdouble dJndxil_leftc=0.0;
gdouble dJpdxipl_leftc=0.0;

gdouble dJdxic_rightc=0.0;
gdouble dJpdxipc_rightc=0.0;
gdouble dJpdphi_rightc=0.0;
gdouble dJdphi_rightc=0.0;

gdouble Bfree=0.0;
gdouble nceq=0.0;
gdouble pceq=0.0;
gdouble Rfree=0.0;

//gdouble nc0_l=0.0;
//gdouble dnc0_l=0.0;
//gdouble pc0_l=0.0;
//gdouble dpc0_l=0.0;

//gdouble nc0_r=0.0;
//gdouble dnc0_r=0.0;
//gdouble pc0_r=0.0;
//gdouble dpc0_r=0.0;

gdouble dJnldxil_l=0.0;
gdouble dJnldxil_c=0.0;
gdouble dJnrdxir_c=0.0;
gdouble dJnrdxir_r=0.0;
gdouble dJpldxipl_l=0.0;
gdouble dJpldxipl_c=0.0;
gdouble dJprdxipr_c=0.0;
gdouble dJprdxipr_r=0.0;

gdouble i0=0.0;
gdouble didv=0.0;  //not a function
gdouble diphic=0.0; //could be a function
gdouble didxic=0.0;
gdouble didxipc=0.0;
gdouble didphir=0.0;
gdouble didxir=0.0;
gdouble didxipr=0.0;
gdouble Nad=0.0;
long double Nion=0.0;

//gdouble dylh_left=0.0;
//gdouble dyrh_left=0.0;
//gdouble dncdphic=0.0;
//gdouble dpcdphic=0.0;
long double dphir_right=0.0;
long double dJndphir_rightr=0.0;
long double dJpdphir_rightr=0.0;
long double dJpdphil_rightc=0.0;
long double dJndphil_rightc=0.0;
long double dJndxir_rightc=0.0;
long double dJpdxipr_rightc=0.0;
long double didphil=0.0;
long double dJdphi_leftc=0.0;
long double dJpdphi_leftc=0.0;
long double dJdxic_leftc=0.0;
long double dJpdxipc_leftc=0.0;
long double didxil=0.0;
long double didxipl=0.0;

long double vl_e=-1.0;
long double vl_h=-1.0;

long double vr_e=-1.0;
long double vr_h=-1.0;


//long double nl0=-1.0;
//long double pl0=-1.0;

//long double p_sh=-1.0;
//long double dp_sh=-1.0;

//long double Jnl_sh=-1.0;
//long double dJnlsh=-1.0;
//long double dJnlsh_c=-1.0;

//long double Jpl_sh=-1.0;
//long double dJplsh=-1.0;
//long double dJplsh_c=-1.0;

struct newton_state *ns=&(in->ns);
struct dimensions *dim=&in->ns.dim;

int contact_left=in->contacts[in->n_contact_y0[z][x]].type;
int contact_right=in->contacts[in->n_contact_y1[z][x]].type;

long double	dJdxic_imag=0.0;
long double dJpdxipc_imag=0.0;

	//if (in->interfaceleft==TRUE)
	//{
	//	ns->phi[z][x][0]=in->Vapplied_y0[z][x];
	//}


		pos=0;
		for (i=0;i<dim->ylen;i++)
		{

			//printf("density=%Le\n",get_dos_ion_density(in,in->imat[z][x][i]));
			//printf("mobility=%Le\n",get_dos_ion_mobility(in,in->imat[z][x][i]));

			if (i==0)
			{

				phil=in->Vapplied_y0[z][x];

				yl=dim->ymesh[0]-(dim->ymesh[1]-dim->ymesh[0]);
				//printf("%Le\n",yl);
				//getchar();
//				Tll=thermal->Ty0;
				Tel=thermal->Ty0;
				Thl=thermal->Ty0;

				Ecl= -in->Xi[z][x][0]-phil;
				Evl= -in->Xi[z][x][0]-phil-in->Eg[z][x][0];
				epl=in->epsilonr[z][x][0]*epsilon0;

				xnl=in->Fi[z][x][0];
				tnl=in->Xi[z][x][0];
				one=xnl+tnl;

				vl_e=in->contacts[in->n_contact_y0[z][x]].ve0;
				nl=get_n_den(in,one,Tel,in->imat[z][x][i]);
				dnl=get_dn_den(in,one,Tel,in->imat[z][x][i]);
				wnl=get_n_w(in,one,Tel,in->imat[z][x][i]);



				munl=in->mun[z][x][0];


				xpl= -in->Fi[z][x][0];
				tpl=(in->Xi[z][x][0]+in->Eg[z][x][0]);
				one=xpl-tpl;

				vl_h=in->contacts[in->n_contact_y0[z][x]].vh0;
				pl=get_p_den(in,one,Thl,in->imat[z][x][i]);
				dpl=get_dp_den(in,one,Thl,in->imat[z][x][i]);
				wpl=get_p_w(in,one,Thl,in->imat[z][x][i]);


				mupl=in->mup[z][x][0];


//				kll=in->kl[i];

			}else
			{
//				Dexl=in->Dex[i-1];
//				exl=in->ex[z][x][i-1];
				phil=ns->phi[z][x][i-1];
				yl=dim->ymesh[i-1];
//				Tll=in->Tl[z][x][i-1];
				Tel=in->Te[z][x][i-1];
				Thl=in->Th[z][x][i-1];
				Ecl=in->Ec[z][x][i-1];
				Evl=in->Ev[z][x][i-1];


				nl=in->n[z][x][i-1];
				dnl=in->dn[z][x][i-1];


				wnl=in->wn[z][x][i-1];
				wpl=in->wp[z][x][i-1];

				pl=in->p[z][x][i-1];
				dpl=in->dp[z][x][i-1];
				munl=in->mun[z][x][i-1];
				mupl=in->mup[z][x][i-1];


				epl=in->epsilonr[z][x][i-1]*epsilon0;


//				kll=in->kl[i-1];
			}

			Ecc=(-in->Xi[z][x][i]-ns->phi[z][x][i]);
			Evc=(-in->Xi[z][x][i]-ns->phi[z][x][i]-in->Eg[z][x][i]);

			if (i==(dim->ylen-1))
			{

//				Dexr=in->Dex[i];
//				exr=0.0;
				//phir=in->V_y1;

				phir=(in->V_y1[z][x]+in->Vapplied_y1[z][x]);



				yr=dim->ymesh[i]+(dim->ymesh[i]-dim->ymesh[i-1]);
//				Tlr=thermal->Ty1;
				Ter=thermal->Ty1;
				Thr=thermal->Ty1;


				Ecr= -in->Xi[z][x][i]-phir;
				Evr= -in->Xi[z][x][i]-phir-in->Eg[z][x][i];


				xnr=(in->V_y1[z][x]+in->Fi[z][x][i]);
				tnr=(in->Xi[z][x][i]);


				one=xnr+tnr;

				vr_e=in->contacts[in->n_contact_y1[z][x]].ve0;
				nr=get_n_den(in,one,Ter,in->imat[z][x][i]);
				dnr=get_dn_den(in,one,Ter,in->imat[z][x][i]);
				wnr=get_n_w(in,one,Ter,in->imat[z][x][i]);



				xpr= -(in->V_y1[z][x]+in->Fi[z][x][i]);
				tpr=(in->Xi[z][x][i]+in->Eg[z][x][i]);

				one=xpr-tpr;

				vr_h=in->contacts[in->n_contact_y1[z][x]].vh0;
				pr=get_p_den(in,one,Thr,in->imat[z][x][i]);
				dpr=get_dp_den(in,one,Thr,in->imat[z][x][i]);
				wpr=get_p_w(in,one,Thr,in->imat[z][x][i]);

				munr=in->mun[z][x][i];
				mupr=in->mup[z][x][i];

				epr=in->epsilonr[z][x][i]*epsilon0;
//				klr=in->kl[i];

			}else
			{

//				Dexr=in->Dex[z][x][i+1];
//				exr=in->ex[z][x][i+1];
				phir=ns->phi[z][x][i+1];
				yr=dim->ymesh[i+1];
//				Tlr=in->Tl[z][x][i+1];
				Ter=in->Te[z][x][i+1];
				Thr=in->Th[z][x][i+1];

				Ecr=in->Ec[z][x][i+1];
				Evr=in->Ev[z][x][i+1];


				nr=in->n[z][x][i+1];
				dnr=in->dn[z][x][i+1];

				wnr=in->wn[z][x][i+1];
				wpr=in->wp[z][x][i+1];

				pr=in->p[z][x][i+1];
				dpr=in->dp[z][x][i+1];
				munr=in->mun[z][x][i+1];
				mupr=in->mup[z][x][i+1];

				epr=in->epsilonr[z][x][i+1]*epsilon0;
//				klr=in->kl[i+1];

			}

			dJdxipc=0.0;
			dJpdxic=0.0;

			epc=in->epsilonr[z][x][i]*epsilon0;


//			exc=in->ex[z][x][i];
//			Dexc=in->Dex[z][x][i];
			yc=dim->ymesh[i];
			dyl=yc-yl;
			dyr=yr-yc;
			ddh=(dyl+dyr)/2.0;
			gdouble dylh=dyl/2.0;
			gdouble dyrh=dyr/2.0;

//			dh=(dyl+dyr);
			phic=ns->phi[z][x][i];
//			Tlc=in->Tl[z][x][i];
//			Tec=in->Te[z][x][i];
//			Thc=in->Th[z][x][i];

				munc=in->mun[z][x][i];
				mupc=in->mup[z][x][i];


				wnc=in->wn[z][x][i];
				wpc=in->wp[z][x][i];

				Dnl=munl*(2.0/3.0)*wnl/Q;
				Dpl=mupl*(2.0/3.0)*wpl/Q;

				Dnc=munc*(2.0/3.0)*wnc/Q;
				Dpc=mupc*(2.0/3.0)*wpc/Q;
				in->Dn[z][x][i]=Dnc;
				in->Dp[z][x][i]=Dnc;


				Dnr=munr*(2.0/3.0)*wnr/Q;
				Dpr=mupr*(2.0/3.0)*wpr/Q;


				Dnl=(Dnl+Dnc)/2.0;
				Dnr=(Dnr+Dnc)/2.0;

				Dpl=(Dpl+Dpc)/2.0;
				Dpr=(Dpr+Dpc)/2.0;

				munl=(munl+munc)/2.0;
				munr=(munr+munc)/2.0;

				mupl=(mupl+mupc)/2.0;
				mupr=(mupr+mupc)/2.0;



				nc=in->n[z][x][i];
				pc=in->p[z][x][i];

				dnc=in->dn[z][x][i];
				dpc=in->dp[z][x][i];
//				dncdphic=in->dndphi[z][x][i];
//				dpcdphic=in->dpdphi[z][x][i];

				Bfree=in->B[z][x][i];
				Nad=in->Nad[z][x][i];
				Nion=in->Nion[z][x][i];

				nceq=in->nfequlib[z][x][i];
				pceq=in->pfequlib[z][x][i];
				Rfree=Bfree*(nc*pc-nceq*pceq);
				in->Rfree[z][x][i]=Rfree;

//			klc=in->kl[i];
			nlast=in->nlast[z][x][i];
			plast=in->plast[z][x][i];

			for (band=0;band<dim->srh_bands;band++)
			{
				in->newton_ntlast[band]=in->ntlast[z][x][i][band];
				in->newton_ptlast[band]=in->ptlast[z][x][i][band];
			}

			dt=in->dt;

//	R=in->R[z][x][i];
	Gn=in->Gn[z][x][i];
	Gp=in->Gp[z][x][i];
	//printf("%d %Le %Le\n",i,in->Gn[z][x][i],in->Gp[z][x][i]);
	e0=(epl+epc)/2.0;
	e1=(epc+epr)/2.0;

//	kl0=(klc+kll)/2.0;
//	kl1=(klr+klc)/2.0;

	dphil= -e0/dyl/ddh;
	dphic= e0/dyl/ddh+e1/dyr/ddh;
	dphir= -e1/dyr/ddh;

	gdouble dphil_d=dphil;
	gdouble dphic_d=dphic;
	gdouble dphir_d=dphir;

	deriv=phil*dphil+phic*dphic+phir*dphir;

	dphidxic=Q*(dnc);
	dphidxipc= -Q*(dpc);


	if (in->ntrapnewton==TRUE)
	{
		for (band=0;band<dim->srh_bands;band++)
		{
			in->newton_dphidntrap[band]=Q*(in->dnt[z][x][i][band]);
			//if ((in->interfaceleft==TRUE)&&(i==0)) in->newton_dphidntrap[band]=0.0;
			//if ((in->interfaceright==TRUE)&&(i==dim->ylen-1)) in->newton_dphidntrap[band]=0.0;
		}
	}

	if (in->ptrapnewton==TRUE)
	{
		for (band=0;band<dim->srh_bands;band++)
		{
			in->newton_dphidptrap[band]= -Q*(in->dpt[z][x][i][band]);
			//if ((in->interfaceleft==TRUE)&&(i==0)) in->newton_dphidptrap[band]=0.0;
			//if ((in->interfaceright==TRUE)&&(i==dim->ylen-1)) in->newton_dphidptrap[band]=0.0;

			//dphidxipc+= -Q*(in->dpt[i]);
		}
	}




//	G=in->G[i];



			xil=Q*2.0*(3.0/2.0)*(Ecc-Ecl)/((wnc+wnl));
			xir=Q*2.0*(3.0/2.0)*(Ecr-Ecc)/((wnr+wnc));

			//gdouble dxil= -Q*2.0*(3.0/2.0)*(Ecc-Ecl)/pow((wnc+wnl),2.0);
			//gdouble dxir= -Q*2.0*(3.0/2.0)*(Ecr-Ecc)/pow((wnr+wnc),2.0);

			xipl=Q*2.0*(3.0/2.0)*(Evc-Evl)/(wpc+wpl);
			xipr=Q*2.0*(3.0/2.0)*(Evr-Evc)/(wpr+wpc);

			dJdxil=0.0;
			dJdxic=0.0;
			dJdxir=0.0;

			dJpdxipl=0.0;
			dJpdxipc=0.0;
			dJpdxipr=0.0;


			dJdphil=0.0;
			dJdphic=0.0;
			dJdphir=0.0;


			dJpdphil=0.0;
			dJpdphic=0.0;
			dJpdphir=0.0;


			Jnl=(Dnl/dyl)*(B(-xil)*nc-B(xil)*nl);
			dJnldxil_l= -(Dnl/dyl)*(B(xil)*dnl);
			dJnldxil_c=(Dnl/dyl)*B(-xil)*dnc;

			gdouble dJnldphi_l= -(munl/dyl)*(dB(-xil)*nc+dB(xil)*nl);
			gdouble dJnldphi_c=(munl/dyl)*(dB(-xil)*nc+dB(xil)*nl);

			Jnr=(Dnr/dyr)*(B(-xir)*nr-B(xir)*nc);
			dJnrdxir_c= -(Dnr/dyr)*(B(xir)*dnc);
			dJnrdxir_r=(Dnr/dyr)*(B(-xir)*dnr);

			gdouble dJnrdphi_c=(munr/dyr)*(-dB(-xir)*nr-dB(xir)*nc);
			gdouble dJnrdphi_r=(munr/dyr)*(dB(-xir)*nr+dB(xir)*nc);

			Jpl=(Dpl/dyl)*(B(-xipl)*pl-B(xipl)*pc);
			dJpldxipl_l=(Dpl/dyl)*(B(-xipl)*dpl);
			dJpldxipl_c= -(Dpl/dyl)*B(xipl)*dpc;

			gdouble dJpldphi_l= -((mupl)/dyl)*(dB(-xipl)*pl+dB(xipl)*pc);
			gdouble dJpldphi_c=((mupl)/dyl)*(dB(-xipl)*pl+dB(xipl)*pc);

			Jpr=(Dpr/dyr)*(B(-xipr)*pc-B(xipr)*pr);
			dJprdxipr_c=(Dpr/dyr)*(B(-xipr)*dpc);
			dJprdxipr_r= -(Dpr/dyr)*(B(xipr)*dpr);

			gdouble dJprdphi_c= -(mupr/dyr)*(dB(-xipr)*pc+dB(xipr)*pr);
			gdouble dJprdphi_r=(mupr/dyr)*(dB(-xipr)*pc+dB(xipr)*pr);


			if (i==0)
			{
				if (contact_left==contact_schottky)
				{
					Jnl=-vl_e*(nc-nl);
					dJnldxil_l= -vl_e*(-dnl);
					dJnldxil_c= -vl_e*(dnc);
					//printf("%Le %Le\n",Jnl,vl_e*(nl0-nl));

					dJnldphi_l= 0.0;//vl_e*(-dnl);
					dJnldphi_c= 0.0;//vl_e*(dnc);

					Jpl=vl_h*(pc-pl);
					dJpldxipl_l= vl_h*(-dpl);
					dJpldxipl_c= vl_h*(dpc);
					//printf("%Le %Le\n",Jnl,vl_e*(nl0-nl));

					dJpldphi_l= 0.0;//-vl_h*(-dpl);//vl_e*(-dnl);
					dJpldphi_c= 0.0;//-vl_h*(dpc);

					dylh=0.0;

				}

				in->Jn_y0[z][x]=Q*Jnl;
				in->Jp_y0[z][x]=Q*Jpl;

			}

			if (i==dim->ylen-1)
			{
				//printf("%d\n",contact_right);
				//getchar();
				if (contact_right==contact_schottky)
				{
					Jnr=-vr_e*(nc-nr);
					dJnrdxir_r= -vr_e*(-dnr);
					dJnrdxir_c= -vr_e*(dnc);
					//printf("%Le %Le\n",Jnl,vr_e*(nl0-nl));

					dJnrdphi_r= 0.0;//vr_e*(-dnl);
					dJnrdphi_c= 0.0;//vr_e*(dnc);

					Jpr=vr_e*(pc-pr);
					dJprdxipr_r= vr_h*(-dpr);
					dJprdxipr_c= vr_h*(dpc);
					//printf("%Le %Le\n",Jnl,vr_e*(nl0-nl));

					dJprdphi_r= 0.0;//vr_e*(-dpl);
					dJprdphi_c= 0.0;//-vr_e*(dpc);

					dyrh=0.0;
				}

				in->Jn_y1[z][x]=Q*Jnr;
				in->Jp_y1[z][x]=Q*Jpr;

			}

			in->Jn[z][x][i]=Q*(Jnl+Jnr)/2.0;
			in->Jp[z][x][i]=Q*(Jpl+Jpr)/2.0;

			dJdxil+= -dJnldxil_l/(dylh+dyrh);
			dJdxic+=(-dJnldxil_c+dJnrdxir_c)/(dylh+dyrh);
			dJdxir+=dJnrdxir_r/(dylh+dyrh);

			dJpdxipl+= -dJpldxipl_l/(dylh+dyrh);
			dJpdxipc+=(-dJpldxipl_c+dJprdxipr_c)/(dylh+dyrh);
			dJpdxipr+=dJprdxipr_r/(dylh+dyrh);


			dJdphil+= -dJnldphi_l/(dylh+dyrh);
			dJdphic+=(-dJnldphi_c+dJnrdphi_c)/(dylh+dyrh);
			dJdphir+=dJnrdphi_r/(dylh+dyrh);


			dJpdphil+= -dJpldphi_l/(dylh+dyrh);
			dJpdphic+=(-dJpldphi_c+dJprdphi_c)/(dylh+dyrh);
			dJpdphir+=dJprdphi_r/(dylh+dyrh);

			if (Bfree!=0.0)
			{
				dJdxic+= -Bfree*(dnc*pc);
				dJdxipc+= -Bfree*(nc*dpc);

				dJpdxipc+=Bfree*(nc*dpc);
				dJpdxic+=Bfree*(dnc*pc);

				//if ((in->interfaceleft==TRUE)&&(i==0))
				//{
				//}else
				//if ((in->interfaceright==TRUE)&&(i==dim->ylen-1))
				//{
				//}else
				//{
					dJdphic+= -Bfree*(dnc*pc);
					dJpdphic+=Bfree*(nc*dpc);
				//}

			}

			if (i==0)
			{
				dJdxic_leftc= dJnldxil_c;
				dJpdxipc_leftc=dJpldxipl_c;
				dJdphi_leftc= dJnldphi_c;
				dJpdphi_leftc=dJpldphi_c;

				dJndphil_leftl=dJnldphi_l;
				dJndphil_leftc=dJnldphi_c;
				dJpdphil_leftl=dJpldphi_l;
				dJpdphil_leftc=dJpldphi_c;

				dphil_left= -e0/dyl/ddh;
				dJndxil_leftc=dJnldxil_c;
				dJpdxipl_leftc=dJpldxipl_c;
				//dylh_left=dylh;
				//dyrh_left=dyrh;

			}

			if (i==dim->ylen-1)
			{
				dJdxic_rightc=dJnrdxir_c;
				dJpdxipc_rightc=dJprdxipr_c;
				dJdphi_rightc= dJnrdphi_c;	//I don't understand the -
				dJpdphi_rightc=dJprdphi_c;

				//new

				dJndphir_rightr=dJnrdphi_r;
				dJndphil_rightc=dJnrdphi_c;
				dJpdphir_rightr=dJprdphi_r;
				dJpdphil_rightc=dJprdphi_c;

				dphir_right= -e1/dyr/ddh;
				dJndxir_rightc=dJnrdxir_c;
				dJpdxipr_rightc=dJprdxipr_c;
			}



			Rtrapn=0.0;
			Rtrapp=0.0;


			in->nrelax[z][x][i]=0.0;
			in->ntrap_to_p[z][x][i]=0.0;
			in->prelax[z][x][i]=0.0;
			in->ptrap_to_n[z][x][i]=0.0;


			if (in->ntrapnewton==TRUE)
			{
				for (band=0;band<dim->srh_bands;band++)
				{
					in->newton_dJdtrapn[band]=0.0;
					in->newton_dJpdtrapn[band]=0.0;
					in->newton_dntrap[band]=nc*in->srh_n_r1[z][x][i][band]-in->srh_n_r2[z][x][i][band]-pc*in->srh_n_r3[z][x][i][band]+in->srh_n_r4[z][x][i][band];
					in->newton_dntrapdntrap[band]=nc*in->dsrh_n_r1[z][x][i][band]-in->dsrh_n_r2[z][x][i][band]-pc*in->dsrh_n_r3[z][x][i][band]+in->dsrh_n_r4[z][x][i][band];
					in->newton_dntrapdn[band]=dnc*in->srh_n_r1[z][x][i][band];
					in->newton_dntrapdp[band]= -dpc*in->srh_n_r3[z][x][i][band];
					Rtrapn+=nc*in->srh_n_r1[z][x][i][band]-in->srh_n_r2[z][x][i][band];
					dJdxic-=dnc*in->srh_n_r1[z][x][i][band];
					in->newton_dJdtrapn[band]-=nc*in->dsrh_n_r1[z][x][i][band]-in->dsrh_n_r2[z][x][i][band];
					Rtrapp+= -(-pc*in->srh_n_r3[z][x][i][band]+in->srh_n_r4[z][x][i][band]);
					dJpdxipc+= -(-dpc*in->srh_n_r3[z][x][i][band]);
					in->newton_dJpdtrapn[band]= -(-pc*in->dsrh_n_r3[z][x][i][band]+in->dsrh_n_r4[z][x][i][band]);


					in->nrelax[z][x][i]+=nc*in->srh_n_r1[z][x][i][band]-in->srh_n_r2[z][x][i][band];
					in->ntrap_to_p[z][x][i]+= -(-pc*in->srh_n_r3[z][x][i][band]+in->srh_n_r4[z][x][i][band]);

					in->nt_r1[z][x][i][band]=nc*in->srh_n_r1[z][x][i][band];
					in->nt_r2[z][x][i][band]=in->srh_n_r2[z][x][i][band];
					in->nt_r3[z][x][i][band]=pc*in->srh_n_r3[z][x][i][band];
					in->nt_r4[z][x][i][band]=in->srh_n_r4[z][x][i][band];

				}
			}

			//band=0;

			if (in->ptrapnewton==TRUE)
			{

				for (band=0;band<dim->srh_bands;band++)
				{
					//dJdtrapn[band]=0.0;
					in->newton_dJpdtrapp[band]=0.0;
					in->newton_dJdtrapp[band]=0.0;
					in->newton_dptrap[band]=pc*in->srh_p_r1[z][x][i][band]-in->srh_p_r2[z][x][i][band]-nc*in->srh_p_r3[z][x][i][band]+in->srh_p_r4[z][x][i][band];
					in->newton_dptrapdptrap[band]=pc*in->dsrh_p_r1[z][x][i][band]-in->dsrh_p_r2[z][x][i][band]-nc*in->dsrh_p_r3[z][x][i][band]+in->dsrh_p_r4[z][x][i][band];
					in->newton_dptrapdp[band]=dpc*in->srh_p_r1[z][x][i][band];
					in->newton_dptrapdn[band]= -dnc*in->srh_p_r3[z][x][i][band];

					Rtrapp+=pc*in->srh_p_r1[z][x][i][band]-in->srh_p_r2[z][x][i][band];
					dJpdxipc+=dpc*in->srh_p_r1[z][x][i][band];
					in->newton_dJpdtrapp[band]+=pc*in->dsrh_p_r1[z][x][i][band]-in->dsrh_p_r2[z][x][i][band];

					Rtrapn+= -(-nc*in->srh_p_r3[z][x][i][band]+in->srh_p_r4[z][x][i][band]);
					dJdxic-= -(-dnc*in->srh_p_r3[z][x][i][band]);
					in->newton_dJdtrapp[band]-= -(-nc*in->dsrh_p_r3[z][x][i][band]+in->dsrh_p_r4[z][x][i][band]);


					in->prelax[z][x][i]+=pc*in->srh_p_r1[z][x][i][band]-in->srh_p_r2[z][x][i][band];
					in->ptrap_to_n[z][x][i]+= -(-nc*in->srh_p_r3[z][x][i][band]+in->srh_p_r4[z][x][i][band]);

					in->pt_r1[z][x][i][band]=pc*in->srh_p_r1[z][x][i][band];
					in->pt_r2[z][x][i][band]=in->srh_p_r2[z][x][i][band];
					in->pt_r3[z][x][i][band]=nc*in->srh_p_r3[z][x][i][band];
					in->pt_r4[z][x][i][band]=in->srh_p_r4[z][x][i][band];

				}

			}

			//band=0;


			in->Rn[z][x][i]=Rtrapn+Rfree;
			in->Rp[z][x][i]=Rtrapp+Rfree;

			in->Rn_srh[z][x][i]=Rtrapn;
			in->Rp_srh[z][x][i]=Rtrapp;
			//Rtrapp=1e24;
			//Rtrapn=1e24;




			if (i!=0)
			{
				mx->Ti[pos]=i;
				mx->Tj[pos]=i-1;
				mx->Tx[pos]=dphil_d;

				pos++;
				//electron
				mx->Ti[pos]=dim->ylen*(1)+i;
				mx->Tj[pos]=dim->ylen*(1)+i-1;
				mx->Tx[pos]=dJdxil;
				pos++;

				mx->Ti[pos]=dim->ylen*(1)+i;
				mx->Tj[pos]=i-1;
				mx->Tx[pos]=dJdphil;
				pos++;

				//hole
				mx->Ti[pos]=dim->ylen*(1+1)+i;
				mx->Tj[pos]=dim->ylen*(1+1)+i-1;
				mx->Tx[pos]=dJpdxipl;
				pos++;

				mx->Ti[pos]=i+dim->ylen*(1+1);
				mx->Tj[pos]=i-1;
				mx->Tx[pos]=dJpdphil;
				pos++;

			}


			//phi
			mx->Ti[pos]=i;
			mx->Tj[pos]=i;
			mx->Tx[pos]=dphic_d;
			pos++;

			mx->Ti[pos]=i;
			mx->Tj[pos]=i+dim->ylen*(1);
			mx->Tx[pos]=dphidxic;
			//strcpy(in->Tdebug[pos],"dphidxic");
			pos++;

			mx->Ti[pos]=i;
			mx->Tj[pos]=i+dim->ylen*(1+1);
			mx->Tx[pos]=dphidxipc;
			//strcpy(in->Tdebug[pos],"dphidxipc");
			pos++;


			//electron

			mx->Ti[pos]=dim->ylen*(1)+i;
			mx->Tj[pos]=dim->ylen*(1)+i;
			mx->Tx[pos]=dJdxic;
			if (mx->complex_matrix==TRUE)
			{
				mx->Txz[pos]=dJdxic_imag;
			}
			//strcpy(in->Tdebug[pos],"dJdxic");
			pos++;


			mx->Ti[pos]=dim->ylen*(1)+i;
			mx->Tj[pos]=dim->ylen*(1+1)+i;
			mx->Tx[pos]=dJdxipc;
			pos++;

			mx->Ti[pos]=dim->ylen*(1)+i;
			mx->Tj[pos]=i;
			mx->Tx[pos]=dJdphic;
			pos++;



			//hole
			mx->Ti[pos]=dim->ylen*(1+1)+i;
			mx->Tj[pos]=dim->ylen*(1+1)+i;
			mx->Tx[pos]=dJpdxipc;
			if (mx->complex_matrix==TRUE)
			{
				mx->Txz[pos]=dJpdxipc_imag;
			}
			pos++;


			mx->Ti[pos]=dim->ylen*(1+1)+i;
			mx->Tj[pos]=dim->ylen*(1)+i;
			mx->Tx[pos]=dJpdxic;
			pos++;

			mx->Ti[pos]=dim->ylen*(1+1)+i;
			mx->Tj[pos]=i;
			mx->Tx[pos]=dJpdphic;
			pos++;


			if (in->ntrapnewton==TRUE)
			{
				for (band=0;band<dim->srh_bands;band++)
				{
					mx->Ti[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tx[pos]=in->newton_dntrapdntrap[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tj[pos]=dim->ylen*1+i;
					mx->Tx[pos]=in->newton_dntrapdn[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tj[pos]=dim->ylen*(1+1)+i;
					mx->Tx[pos]=in->newton_dntrapdp[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tx[pos]=in->newton_dJdtrapn[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tx[pos]=in->newton_dJpdtrapn[band];
					pos++;

					mx->Ti[pos]=i;
					mx->Tj[pos]=dim->ylen*(1+1+1+band)+i;
					mx->Tx[pos]=in->newton_dphidntrap[band];
					pos++;

				}


			}

			if (in->ptrapnewton==TRUE)
			{
				for (band=0;band<dim->srh_bands;band++)
				{
					mx->Ti[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tx[pos]=in->newton_dptrapdptrap[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tj[pos]=dim->ylen*(1+1)+i;
					mx->Tx[pos]=in->newton_dptrapdp[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tj[pos]=dim->ylen*(1)+i;
					mx->Tx[pos]=in->newton_dptrapdn[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1+1)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tx[pos]=in->newton_dJpdtrapp[band];
					pos++;

					mx->Ti[pos]=dim->ylen*(1)+i;
					mx->Tj[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tx[pos]=in->newton_dJdtrapp[band];
					pos++;

					mx->Ti[pos]=i;
					mx->Tj[pos]=dim->ylen*(1+1+1+dim->srh_bands+band)+i;
					mx->Tx[pos]=in->newton_dphidptrap[band];
					pos++;
				}

			}

			if (i!=(dim->ylen-1))
			{


				//phi
				mx->Ti[pos]=i;
				mx->Tj[pos]=i+1;
				mx->Tx[pos]=dphir_d;
				pos++;


				//electron
				mx->Ti[pos]=dim->ylen*(1)+i;
				mx->Tj[pos]=dim->ylen*(1)+i+1;
				mx->Tx[pos]=dJdxir;
				pos++;

				mx->Ti[pos]=i+dim->ylen*(1);
				mx->Tj[pos]=i+1;
				mx->Tx[pos]=dJdphir;
				pos++;

				//hole
				mx->Ti[pos]=dim->ylen*(1+1)+i;
				mx->Tj[pos]=dim->ylen*(1+1)+i+1;
				mx->Tx[pos]=dJpdxipr;
				pos++;

				mx->Ti[pos]=i+dim->ylen*(1+1);
				mx->Tj[pos]=i+1;
				mx->Tx[pos]=dJpdphir;
				pos++;


			}

			//Possion
			gdouble build=0.0;

			build= -(deriv);

			build+= -(-(pc-nc+Nad+Nion)*Q);

			for (band=0;band<dim->srh_bands;band++)
			{
				build+= -(-Q*(in->pt[z][x][i][band]-in->nt[z][x][i][band]));
			}

			//build+= -(-Q*in->Nad[i]);

			mx->b[i]=build;


			//Electron
			build=0.0;
			build= -((Jnr-Jnl)/(dylh+dyrh)-Rtrapn-Rfree);


			//getchar();
			build-=Gn;
			mx->b[dim->ylen*(1)+i]=build;

			if (mx->complex_matrix==TRUE)
			{
				build=0.0;
				build-=in->omega*Gn*0.05;
				mx->bz[dim->ylen*(1)+i] = build;
			}

			//hole
			build=0.0;
			build= -((Jpr-Jpl)/(dylh+dyrh)+Rtrapp+Rfree);

			build-= -Gp;


			mx->b[dim->ylen*(1+1)+i]=build;

			if (mx->complex_matrix==TRUE)
			{
				build=0.0;
				build-=-in->omega*Gp*0.05;
				mx->bz[dim->ylen*(1+1)+i] = build;
			}

			if (in->ntrapnewton==TRUE)
			{
				for (band=0;band<dim->srh_bands;band++)
				{
					mx->b[dim->ylen*(1+1+1+band)+i]= -(in->newton_dntrap[band]);
				}
			}

			if (in->ptrapnewton==TRUE)
			{
				for (band=0;band<dim->srh_bands;band++)
				{
					mx->b[dim->ylen*(1+1+1+dim->srh_bands+band)+i]= -(in->newton_dptrap[band]);
				}

			}

		}

if (pos>mx->nz)
{
	ewe(sim,"Error %d %d %d\n",pos,mx->nz,in->kl_in_newton);
}

}

gdouble get_cur_error(struct simulation *sim,struct device *in)
{
int i;
gdouble phi=0.0;
gdouble n=0.0;
gdouble p=0.0;
gdouble x=0.0;
gdouble te=0.0;
gdouble th=0.0;
gdouble tl=0.0;
gdouble ttn=0.0;
gdouble ttp=0.0;
gdouble i0=0.0;
int band=0;
struct dimensions *dim=&in->ns.dim;
struct matrix *mx = &(in->mx);

for (i=0;i<dim->ylen;i++)
{
		//if ((in->interfaceleft==TRUE)&&(i==0))
		//{
		//}else
		//if ((in->interfaceright==TRUE)&&(i==dim->ylen-1))
		//{
		//}else
		//{
			phi+=gfabs(mx->b[i]);
		//}

		n+=gfabs(mx->b[dim->ylen*(1)+i]);
		p+=+gfabs(mx->b[dim->ylen*(1+1)+i]);

		if (in->ntrapnewton==TRUE)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				ttn+=gfabs(mx->b[dim->ylen*(1+1+1+band)+i]);
			}
		}

		if (in->ptrapnewton==TRUE)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				ttp+=gfabs(mx->b[dim->ylen*(1+1+1+dim->srh_bands+band)+i]);
			}
		}


}


gdouble tot=phi+n+p+x+te+th+tl+ttn+ttp+i0;
if (isnan( tot))
{
	printf_log(sim,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",phi,n,p,x,te,th,tl,ttn,ttp);
	//dump_matrix(in);
	ewe(sim,"nan detected in newton solver\n");
}

return tot;
}

void solver_cal_memory(struct device *in,int *ret_N,int *ret_M)
{
int i=0;
int N=0;
int M=0;
struct dimensions *dim=&in->ns.dim;
//struct matrix *mx=&(in->mx);


//in->contact_shift_l=0;


N=dim->ylen*3-2;	//Possion main

N+=dim->ylen*3-2;	//Je main
N+=dim->ylen*3-2;	//Jh main

N+=dim->ylen*3-2;	//dJe/phi
N+=dim->ylen*3-2;	//dJh/phi

N+=dim->ylen;		//dphi/dn
N+=dim->ylen;		//dphi/dh

N+=dim->ylen;		//dJndp
N+=dim->ylen;		//dJpdn

M=dim->ylen;			//Pos
M+=dim->ylen;		//Je
M+=dim->ylen;		//Jh

if (in->ntrapnewton==TRUE)
{
	for (i=0;i<dim->srh_bands;i++)
	{
		N+=dim->ylen;		//dntrapdn
		N+=dim->ylen;		//dntrapdntrap
		N+=dim->ylen;		//dntrapdp
		N+=dim->ylen;		//dJndtrapn
		N+=dim->ylen;		//dJpdtrapn
		N+=dim->ylen;		//dphidntrap

		M+=dim->ylen;		//nt
	}

}

if (in->ptrapnewton==TRUE)
{
	for (i=0;i<dim->srh_bands;i++)
	{
		N+=dim->ylen;		//dptrapdp
		N+=dim->ylen;		//dptrapdptrap
		N+=dim->ylen;		//dptrapdn
		N+=dim->ylen;		//dJpdtrapp
		N+=dim->ylen;		//dJdtrapp
		N+=dim->ylen;		//dphidptrap

		M+=dim->ylen;		//pt
	}
}


*ret_N=N;
*ret_M=M;
}

void dllinternal_solver_realloc(struct simulation *sim,struct device *in, int dim_)
{
int N=0;
int M=0;
struct dimensions *dim=&in->ns.dim;
struct matrix *mx=&(in->mx);
gdouble *dtemp=NULL;



solver_cal_memory(in,&N,&M);

//printf("realloc %d %d\n",mx->nz,N);
//getchar();

int alloc=FALSE;
if ((mx->nz==0)||(mx->M==0))
{
	mx->nz=N;
	mx->M=M;
	alloc=TRUE;
}else
if ((N!=mx->nz)||(M!=mx->M))
{
	mx->nz=N;
	mx->M=M;
	alloc=TRUE;
}


if (alloc==TRUE)
{

	matrix_realloc(sim,mx);

	if (dim->srh_bands>0)
	{
		dtemp=realloc(in->newton_dntrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dntrap=dtemp;
		}

		dtemp=realloc(in->newton_dntrapdntrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dntrapdntrap=dtemp;
		}

		dtemp=realloc(in->newton_dntrapdn,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dntrapdn=dtemp;
		}

		dtemp=realloc(in->newton_dntrapdp,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dntrapdp=dtemp;
		}

		dtemp=realloc(in->newton_dJdtrapn,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dJdtrapn=dtemp;
		}

		dtemp=realloc(in->newton_dJpdtrapn,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dJpdtrapn=dtemp;
		}

		dtemp=realloc(in->newton_dphidntrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dphidntrap=dtemp;
		}


		dtemp=realloc(in->newton_ntlast,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_ntlast=dtemp;
		}



		dtemp=realloc(in->newton_dptrapdp,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dptrapdp=dtemp;
		}

		dtemp=realloc(in->newton_dptrapdptrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dptrapdptrap=dtemp;
		}

		dtemp=realloc(in->newton_dptrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dptrap=dtemp;
		}

		dtemp=realloc(in->newton_dptrapdn,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dptrapdn=dtemp;
		}

		dtemp=realloc(in->newton_dJpdtrapp,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dJpdtrapp=dtemp;
		}


		dtemp=realloc(in->newton_dJdtrapp,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dJdtrapp=dtemp;
		}

		dtemp=realloc(in->newton_dphidptrap,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_dphidptrap=dtemp;
		}

		dtemp=realloc(in->newton_ptlast,dim->srh_bands*sizeof(gdouble));
		if (dtemp==NULL)
		{
			ewe(sim,"memory error\n");
		}else
		{
			in->newton_ptlast=dtemp;
		}

	}

}

}

void dllinternal_solver_free_memory(struct simulation *sim,struct device *in)
{
//int i=0;
struct dimensions *dim=&in->ns.dim;
if (dim->srh_bands>0)
{
	free(in->newton_dntrap);
	free(in->newton_dntrapdntrap);
	free(in->newton_dntrapdn);
	free(in->newton_dntrapdp);
	free(in->newton_dJdtrapn);
	free(in->newton_dJpdtrapn);
	free(in->newton_dphidntrap);
	free(in->newton_dptrapdp);
	free(in->newton_dptrapdptrap);
	free(in->newton_dptrap);
	free(in->newton_dptrapdn);
	free(in->newton_dJpdtrapp);
	free(in->newton_dJdtrapp);
	free(in->newton_dphidptrap);

	free(in->newton_ntlast);
	free(in->newton_ptlast);

}


//for (i=0;i<mx->nz;i++)
//{
//	free(in->Tdebug[i]);
//}
//free(in->Tdebug);

in->newton_dntrapdntrap=NULL;
in->newton_dntrap=NULL;
in->newton_dntrapdn=NULL;
in->newton_dntrapdp=NULL;
in->newton_dJdtrapn=NULL;
in->newton_dJpdtrapn=NULL;
in->newton_dphidntrap=NULL;
in->newton_dptrapdp=NULL;
in->newton_dptrapdptrap=NULL;
in->newton_dptrap=NULL;
in->newton_dptrapdn=NULL;
in->newton_dJpdtrapp=NULL;
in->newton_dJdtrapp=NULL;
in->newton_dphidptrap=NULL;


in->newton_ntlast=NULL;
in->newton_ptlast=NULL;

matrix_free(sim,&(in->mx));
//in->Tdebug=NULL;

}

int dllinternal_solve_cur(struct simulation *sim,struct device *in, int z, int x)
{
gdouble error=0.0;
int ittr=0;
char temp[PATH_MAX];
struct dimensions *dim=&in->ns.dim;
struct newton_state *ns=&(in->ns);
struct matrix *mx=&(in->mx);

//if (get_dump_status(sim,dump_print_newtonerror)==TRUE)
//{
//	printf_log(sim,"Solve cur\n");
//}


#ifdef only
only_update_thermal=FALSE;
#endif
//in->enable_back=FALSE;
int stop=FALSE;
int thermalrun=0;
gdouble check[10];
int cpos=0;
mx->use_cache=in->cache.enabled;

matrix_cache_reset(sim,mx);
	do
	{
		fill_matrix(sim,in,z,x);

		if (in->newton_only_fill_matrix==TRUE)
		{
			return 0;
		}


		//dump_for_plot(in);
		//plot_now(in,"plot");
			//solver_dump_matrix(in->M,mx->nz,in->Ti,in->Tj, in->Tx,mx->b);
		//getchar();

		if (in->stop==TRUE)
		{
			break;
		}
		//matrix_solve(sim,mx);

		if (matrix_solve(sim,mx)==0)
		{
			//printf("found in cache %s\n",mx->cache_file_path);
			//getchar();
			newton_state_load(sim,ns,mx->cache_file_path);
			newton_state_update_device(sim,in, ns);
			//getchar();
			return 0;
		}/*else
		{
			printf("state not foudn %s\n",mx->cache_file_path);
			getchar();
		}*/

		update_solver_vars(sim,in,z,x,TRUE);

		//solver_dump_matrix(in->M,mx->nz,in->Ti,in->Tj, in->Tx,mx->b);
		//getchar();


		error=get_cur_error(sim,in);

		//thermalrun++;
		if (thermalrun==40) thermalrun=0;
		//update(in);
//getchar();

		if (get_dump_status(sim,dump_print_newtonerror)==TRUE)
		{
			printf_log(sim,"%d Cur error = %Le %Le I=%Le",ittr,error,contact_get_active_contact_voltage(sim,in),get_I(in));

			printf_log(sim,"\n");
			//printf_log(sim,"%Le\n",get_equiv_V(sim,in));
		}

		ns->last_error=error;
		ns->last_ittr=ittr;
		ittr++;

		if (get_dump_status(sim,dump_write_converge)==TRUE)
		{
			sim->converge=fopena(get_output_path(sim),"converge.dat","a");
			fprintf(sim->converge,"%Le\n",error);
			fclose(sim->converge);
		}

		stop=TRUE;

		if (ittr<in->max_electrical_itt)
		{
			if (error>in->min_cur_error)
			{
				stop=FALSE;
			}
		}

		if (ittr<in->newton_min_itt)
		{
			stop=FALSE;
		}


		if (in->newton_clever_exit==TRUE)
		{
			check[cpos]=error;
			cpos++;

			if (cpos>10)
			{
				cpos=0;
			}

			if (ittr>=in->newton_min_itt)
			{
					if ((check[0]<error)||(check[1]<error))
					{
						stop=TRUE;
					}
			}
		}

		if ((ittr<2)&&(error<in->min_cur_error))
		{
			in->dd_conv=TRUE;
		}else
		{
			in->dd_conv=FALSE;
		}

	}while(stop==FALSE);

in->newton_last_ittr=ittr;

if (error>1e-3)
{
	printf_log(sim,"warning: The solver has not converged very well.\n");
}


if (get_dump_status(sim,dump_newton)==TRUE)
{
	join_path(2,temp,get_output_path(sim),"solver");
	dump_1d_slice(sim,in,temp);
}
//plot_now(sim,in,"plot");
//getchar();
in->odes+=mx->M;

if (mx->use_cache==TRUE)
{
	newton_state_save(sim,mx->cache_file_path,ns);
}
//getchar();

return 0;
}


