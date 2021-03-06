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

/** @file gendosfdgause.c
	@brief Generates the DoS files.
*/

#include <stdio.h>
#include <stdlib.h>

	#include <zlib.h>

#include <code_ctrl.h>
#include <sim.h>
#include <inp.h>
#include <util.h>
#include <dat_file.h>
#include <hard_limit.h>
#include <epitaxy.h>
#include <lang.h>
#include <dos.h>

static int unused __attribute__((unused));

#include "dos_an.h"

#include "server.h"
#include "dump.h"
#include "log.h"
#include "checksum.h"
#include "cal_path.h"

#ifdef dos_debug
#define test_dist


gdouble pick[1000];
gdouble pick_x[1000];
gdouble pick_start= -7;
gdouble pick_stop=0;
gdouble pick_dx;
#endif

static struct dosconfig confige[10];
static struct dosconfig configh[10];

void dos_double_res(struct simulation *sim)
{
/*gdouble number;
FILE *in=fopen("srhbandp.inp","r");
FILE *out=fopen("srhbandp.new","w");
do
{
unused=fscanf(in,"%Le",&number);
fprintf(out,"%Le\n",number);
fprintf(out,"%Le\n",number);
}while(!feof(in));
fclose(in);
fclose(out);
*/
int i;
for (i=1;i<=80;i++)
{
printf_log(sim,"p srhbandp.inp %d 1\n",i);
}
}

void dump_qe(struct simulation *sim)
{
struct math_xy n;
struct math_xy p;

inter_load(sim,&n,"dosoutn.dat");
inter_load(sim,&p,"dosoutp.dat");
inter_swap(&p);


gdouble start=n.x[0];
gdouble stop=n.x[n.len-1];
gdouble step=5e-3;
gdouble pos=0.0;
gdouble hf=0.0;
gdouble hf_step=1e-2;
gdouble re=0.0;
gdouble rh=0.0;
gdouble rhored=0.0;
gdouble rhored_tot=0.0;
//gdouble Ef=Ec-0.5;
pos=start;
gdouble max=0.0;
FILE *qe=fopen("qe.dat","w");
hf_step=1e-1;
do
{
pos=start;
rhored_tot=0.0;
	do
	{
	re=inter_get(&n,pos+hf);
	rh=inter_get(&p,pos);
	rhored=re*rh;

	rhored_tot+=rhored;
	pos+=step;
	}while(pos<stop);

if (rhored_tot>max) max=rhored_tot;
//getchar();
hf+=hf_step;
//printf_log(sim,"%Le\n",hf);
}while(hf<2.0);

hf=0;
hf_step=1e-2;

do
{
pos=start;
rhored_tot=0.0;
	do
	{
	re=inter_get(&n,pos+hf);
	rh=inter_get(&p,pos);
	rhored=re*rh;
	rhored_tot+=rhored;

	pos+=step;
	}while(pos<stop);
//getchar();
fprintf(qe,"%Le %Le\n",hf,rhored_tot/max);
hf+=hf_step;
//printf_log(sim,"%Le\n",hf);
}while(hf<2.0);

fclose(qe);
inter_free(&n);
inter_free(&p);
}

void pick_init()
{
#ifdef dos_debug
int i;
pick_dx=(pick_stop-pick_start)/1000.0;
gdouble pos=pick_start;
for (i=0;i<1000;i++)
{
pick[i]=0.0;
pick_x[i]=pos;
pos+=pick_dx;
}
#endif
}

void pick_add(gdouble x,gdouble value)
{
#ifdef dos_debug
int pos;
for (pos=0;pos<1000;pos++)
{
if (pick_x[pos]>x) break;
}
if (pos<0.0) return;
if (pos>999) return;
if (value>pick[pos]) pick[pos]=value;
#endif
}

void pick_dump()
{
#ifdef dos_debug
int i;
FILE *out;
out=fopen("gaus.dat","w");
gdouble pos=pick_start;
for (i=0;i<1000;i++)
{
fprintf(out,"%Le %Le\n",pick_x[i],pick[i]);
pos+=pick_dx;
}
fclose(out);
#endif
}

void gen_do(struct simulation *sim,struct dosconfig *in,struct dosconfig *in2,char * outfile,int electrons,int mat)
{
char name[200];
char temp[1000];
gdouble tstart=0.0;
gdouble tstop=0.0;
gdouble tsteps=0.0;
gdouble dt=0.0;
gdouble xpos=0.0;
gdouble tpos=0.0;
int t=0;
int x=0;
int band;
gdouble *srh_r1=NULL;
gdouble *srh_r2=NULL;
gdouble *srh_r3=NULL;
gdouble *srh_r4=NULL;
gdouble *srh_n=NULL;
gdouble *srh_den=NULL;
gdouble *srh_dE_sum=NULL;
gdouble *srh_read=NULL;

gdouble *srh_x=NULL;
gdouble *srh_mid=NULL;
gdouble *band_E_mesh=NULL;
int *band_i=NULL;

int e=0;
gdouble E=0.0;
gdouble dE=0.0;
gdouble rho=0.0;
gdouble sum=0.0;
gdouble f=0.0;
gdouble last_n0=0;
gdouble *xmesh=NULL;

long double srh_pos=0.0;
long double srh_delta=0.0;

int i;
int band_pos=0;
int cur_band=0;
int points_per_band=0;

long double E_free_start=0.0;
long double E_free_stop=2.0;
int E_free_points=1000.0;
long double dE_free=(E_free_stop-E_free_start)/(long double)E_free_points;
long double sum_l=0.0;
long double sum_r=0.0;
long double pos=0.0;
struct dos_an_data my_dos_an;

//#define dos_test_stats
#ifdef dos_test_stats
FILE *freetest;
if (electrons==TRUE)
{
	sprintf(name,"freetestn.dat");
	remove_file(sim,name);
}else
{
	sprintf(name,"freetestp.dat");
	remove_file(sim,name);
}
#endif

if (in->srh_bands!=0)
{
	srh_r1=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_r2=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_r3=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_r4=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_n=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_den=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_dE_sum=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_read=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_x=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
	srh_mid=(gdouble *)malloc(sizeof(gdouble)*in->srh_bands);
}

if (in->Esteps!=0)
{
	band_E_mesh=(gdouble *)malloc(sizeof(gdouble)*in->Esteps);
	band_i=(int *)malloc(sizeof(int)*in->Esteps);
}


tstart=in->Tstart;
tstop=in->Tstop;
tsteps=in->Tsteps;
dt=(tstop-tstart)/tsteps;

gdouble dxr=(in->nstop-in->nstart)/((gdouble)in->npoints);

xpos=in->nstart;
tpos=tstart;
t=0;
x=0;

int buf_len=0;
buf_len+=19;
buf_len+=in->npoints;		//mesh
buf_len+=tsteps;		//mesh
buf_len+=in->srh_bands;	//mesh
buf_len+=in->srh_bands;	//holds the density
buf_len+=tsteps*in->npoints*2; //data
buf_len+=tsteps*in->npoints*5*in->srh_bands; //data

int buf_pos=0;
gdouble *buf=(gdouble*)malloc(buf_len*sizeof(gdouble));
buf[buf_pos++]=(gdouble)in->npoints;
buf[buf_pos++]=(gdouble)tsteps;
buf[buf_pos++]=(gdouble)in->srh_bands;
buf[buf_pos++]=(gdouble)in->epsilonr;
buf[buf_pos++]=(gdouble)in->doping_start;
buf[buf_pos++]=(gdouble)in->doping_stop;
buf[buf_pos++]=(gdouble)in->mu;
buf[buf_pos++]=(gdouble)in->ion_density;
buf[buf_pos++]=(gdouble)in->ion_mobility;
buf[buf_pos++]=(gdouble)in->srh_vth;
buf[buf_pos++]=(gdouble)in->srh_sigman;
buf[buf_pos++]=(gdouble)in->srh_sigmap;
buf[buf_pos++]=(gdouble)in->Nc;
buf[buf_pos++]=(gdouble)in->Nv;
buf[buf_pos++]=(gdouble)in->Nt;
buf[buf_pos++]=(gdouble)in->Eg;
buf[buf_pos++]=(gdouble)in->Xi;
buf[buf_pos++]=(gdouble)in->B;
buf[buf_pos++]=(gdouble)in->dos_free_carrier_stats;

xmesh=(gdouble *)malloc(sizeof(gdouble)*in->npoints);

for (x=0;x<in->npoints;x++)
{
	buf[buf_pos++]=xpos;
	xmesh[x]=xpos;
	xpos+=dxr;
}

for (t=0;t<tsteps;t++)
{
	buf[buf_pos++]=tpos;
	tpos+=dt;
}

//Allocate where the traps are going to be in energy space
srh_delta=fabs(in->srh_stop-in->srh_start)/(long double)(in->srh_bands);
srh_pos=in->srh_start;

for (band=0;band<in->srh_bands;band++)
{
	srh_pos+=srh_delta/2.0;
	srh_mid[band]=srh_pos;
	buf[buf_pos++]=srh_pos;

	srh_pos+=srh_delta/2.0;
	srh_x[band]=srh_pos;

}


//Build the energy mesh for integration
if (in->srh_bands>0)
{
	points_per_band=in->Esteps/in->srh_bands;
}

pos=in->srh_start;
dE=fabs(in->srh_stop-in->srh_start)/((gdouble)in->Esteps);

//FILE *test4=fopen("test4.dat","w");
for (i=0;i<in->Esteps;i++)
{
	band_E_mesh[i]=pos;
	band_i[i]=cur_band;
	//fprintf(test4,"%Le %d\n",band_E_mesh[i],band_i[i]);
	//getchar();
	pos+=dE;
	band_pos++;
	if (band_pos>=points_per_band)
	{
		band_pos=0;
		cur_band++;
	}



}
//fclose(test4);
//getchar();

tpos=tstart;

if (in->dostype==dos_read)
{
FILE *dosread;
	if (electrons==TRUE)
	{
		sprintf(name,"%s_srhbandn.dat",confige[mat].dos_name);
		dosread=fopen(name,"r");
	}else
	{
		sprintf(name,"%s_srhbandp.dat",configh[mat].dos_name);
		dosread=fopen(name,"r");
	}

	if (dosread==NULL)
	{
		ewe(sim,"can not open srhband file\n");
	}

	for (band=0;band<in->srh_bands;band++)
	{
		unused=fscanf(dosread,"%Le",&(srh_read[band]));
		printf_log(sim,"%Le\n",srh_read[band]);
		srh_read[band]=fabs(srh_read[band]);
	}


	fclose(dosread);
}else
if (in->dostype==dos_an)
{
printf_log(sim,"---------------------Generating---------------------\n");
	if (electrons==TRUE)
	{
		sprintf(name,"%s.inp",confige[mat].analytical_dos_file_name);
		dos_an_load(sim,&my_dos_an,name);
	}else
	{
		sprintf(name,"%s.inp",configh[mat].analytical_dos_file_name);
		dos_an_load(sim,&my_dos_an,name);
	}

	for (band=0;band<in->srh_bands;band++)
	{
		srh_read[band]=fabs(dos_an_get_value(sim,&my_dos_an,srh_mid[band]));
		printf_log(sim,"%d %Le\n",band,srh_read[band]);
	}

}

#ifdef test_dist
FILE *munfile=fopen("munfile.dat","w");
if (munfile==NULL)
{
	ewe("problem\n");
}


FILE *plotbands;
plotbands=fopen("plotbandsn.dat","w");

FILE *plotbandsfree;
plotbandsfree=fopen("plotbandsfreen.dat","w");




FILE *rod=fopen("gau_test_n.dat","w");

#endif


#ifdef test_dist
int first=FALSE;
#endif

struct dat_file dos_out;
buffer_init(&dos_out);

if (get_dump_status(sim,dump_write_out_band_structure)==TRUE)
{
	if (electrons==TRUE)
	{
		buffer_malloc(&dos_out);
		dos_out.y_mul=1.0;
		dos_out.x_mul=1.0;
		strcpy(dos_out.title,"Discretized density of states");
		strcpy(dos_out.type,"xy");
		strcpy(dos_out.x_label,"Energy");
		strcpy(dos_out.y_label,"Density of states");
		strcpy(dos_out.x_units,"eV");
		strcpy(dos_out.y_units,"m^{-3}eV^{-1}");
		dos_out.logscale_x=0;
		dos_out.logscale_y=0;
		buffer_add_info(sim,&dos_out);
	}else
	{
		buffer_malloc(&dos_out);
		dos_out.y_mul=1.0;
		dos_out.x_mul=1.0;
		strcpy(dos_out.title,"Discretized density of states");
		strcpy(dos_out.type,"xy");
		strcpy(dos_out.x_label,"Energy");
		strcpy(dos_out.y_label,"Density of states");
		strcpy(dos_out.x_units,"eV");
		strcpy(dos_out.y_units,"m^{-3}eV^{-1}");
		dos_out.logscale_x=0;
		dos_out.logscale_y=0;
		buffer_add_info(sim,&dos_out);
	}
}

int srh_band=0;
gdouble srh_E=0.0;
gdouble srh_f=0.0;
for (t=0;t<tsteps;t++)
{

printf_log(sim,"%d/%d\n",t,(int)tsteps);

	for (x=0;x<in->npoints;x++)
	{
		xpos=xmesh[x];
		E=0.0;

		for (srh_band=0;srh_band<(in->srh_bands);srh_band++)
		{
			srh_r1[srh_band]=0.0;
			srh_r2[srh_band]=0.0;
			srh_r3[srh_band]=0.0;
			srh_r4[srh_band]=0.0;
			srh_n[srh_band]=0.0;
			srh_den[srh_band]=0.0;
			srh_dE_sum[srh_band]=0.0;
		}

		srh_band=0;

		for (e=0;e<in->Esteps;e++)
		{
			E=band_E_mesh[e];
			srh_band=band_i[e];
			srh_E=srh_mid[srh_band];

			f=1.0/(1.0+exp((E-xpos)*Q/(kb*tpos)));

			srh_f=1.0/(1.0+exp((srh_E-xpos)*Q/(kb*tpos)));

			if (in->dostype==dos_exp)
			{
				rho=in->Nt*exp((E)/(in->Et));
			}else
			if ((in->dostype==dos_read)||(confige[mat].dostype==dos_an))
			{
				rho=srh_read[srh_band];
			}else
			{
				ewe(sim,"I don't understand this DoS type.\n");
			}


			if (x==0)
			{
				pick_add(E-in->Xi,rho);
				if (get_dump_status(sim,dump_write_out_band_structure)==TRUE)
				{
					if (E>in->srh_start)
					{
						if (electrons==TRUE)
						{
							sprintf(temp,"%Le %Le\n",E-in->Xi,rho);
							buffer_add_string(&dos_out,temp);
						}else
						{
							sprintf(temp,"%Le %Le\n",-E-in->Xi-in->Eg,rho);
							buffer_add_string(&dos_out,temp);
						}

					}
				}
			}


			if (electrons==TRUE)
			{
				srh_r1[srh_band]+=in->srh_vth*in->srh_sigman*rho*(1.0-srh_f)*dE;
				srh_r2[srh_band]+=in->srh_vth*in->srh_sigman*in->Nc*exp((srh_E*Q)/(tpos*kb))*rho*srh_f*dE;//in->srh_vth*srh_sigman*rho*(1.0-srh_f)*dE;//
				srh_r3[srh_band]+=in->srh_vth*in->srh_sigmap*rho*srh_f*dE;
				srh_r4[srh_band]+=in->srh_vth*in->srh_sigmap*in->Nv*exp((-in->Eg*Q-srh_E*Q)/(tpos*kb))*rho*(1.0-srh_f)*dE;//in->srh_vth*srh_sigmap*rho*srh_f*dE;//
				//if (srh_r1[srh_band]<1e-3)
				//{
				//	printf("%Lf %Le %Le %Le %Le\n",srh_mid[srh_band],srh_r1[srh_band],srh_r2[srh_band],srh_r3[srh_band],srh_r4[srh_band]);
				//}
				srh_n[srh_band]+=rho*srh_f*dE;
				srh_den[srh_band]+=rho*dE;
				srh_dE_sum[srh_band]+=dE;
			}else
			{
				srh_r1[srh_band]+=in->srh_vth*in->srh_sigmap*rho*(1.0-srh_f)*dE;
				srh_r2[srh_band]+=in->srh_vth*in->srh_sigmap*in->Nv*exp((srh_E*Q)/(tpos*kb))*rho*srh_f*dE;//in->srh_vth*srh_sigman*rho*(1.0-srh_f)*dE;//
				srh_r3[srh_band]+=in->srh_vth*in->srh_sigman*rho*srh_f*dE;
				srh_r4[srh_band]+=in->srh_vth*in->srh_sigman*in->Nc*exp((-in->Eg*Q-(srh_E*Q))/(tpos*kb))*rho*(1.0-srh_f)*dE;//in->srh_vth*srh_sigmap*rho*srh_f*dE;//
				srh_n[srh_band]+=rho*srh_f*dE;
				srh_den[srh_band]+=rho*dE;
				srh_dE_sum[srh_band]+=dE;
			}

		}



		if (x==0)
		{
			for (band=0;band<in->srh_bands;band++)
			{
				buf[buf_pos++]=srh_den[band];
			}

			if (get_dump_status(sim,dump_band_structure)==TRUE)
			{
				FILE *bandsdump;
				if (electrons==TRUE)
				{
					bandsdump=fopen("lumo_out.dat","w");
					for (band=0;band<in->srh_bands;band++)
					{
						fprintf(bandsdump,"%Le\n",srh_den[band]/srh_dE_sum[band]);
					}
					fclose(bandsdump);
				}else
				{
					bandsdump=fopen("homo_out.dat","w");
					for (band=0;band<in->srh_bands;band++)
					{
						fprintf(bandsdump,"%Le\n",srh_den[band]/srh_dE_sum[band]);
					}
					fclose(bandsdump);
				}
			}
		}

		///////////Free stuff
		sum_l=0.0;
		sum_r=0.0;
		long double Nc=2.0*powl((in->me*m0*kb*tpos)/(2.0*PI*hbar*hbar),3.0/2.0);
		long double Nv=2.0*powl((in->mh*m0*kb*tpos)/(2.0*PI*hbar*hbar),3.0/2.0);

		if (in->dos_free_carrier_stats==fd_look_up_table)
		{
			E=E_free_start;
			sum=0.0;

			gdouble w0=0.0;
			for (e=0;e<E_free_points;e++)
			{

				if (electrons==TRUE)
				{
					rho=(sqrtl(E*Q)/(2.0*PI*PI))*pow((2.0*in->me*m0)/(hbar*hbar),3.0/2.0);
				}else
				{
					rho=(sqrtl(E*Q)/(2.0*PI*PI))*pow((2.0*in->mh*m0)/(hbar*hbar),3.0/2.0);
				}

				f=1.0/(1.0+expl((E-xpos)*Q/(kb*tpos)));
				sum+=rho*Q*f*dE_free;

				f=1.0/(1.0+expl((E-xpos-dxr/2.0)*Q/(kb*tpos)));
				sum_l+=rho*Q*f*dE_free;

				f=1.0/(1.0+expl((E-xpos+dxr/2.0)*Q/(kb*tpos)));
				sum_r+=rho*Q*f*dE_free;
				E+=dE_free;
				//if (electrons==FALSE)
				//{
				//	printf("%Le %Le\n",rho,in->mh);
				//	getchar();
				//}

			}

			//printf("here %Le %Le %Le \n",sum,in->Nc*exp((xpos*Q)/(kb*tpos)),xpos);//
			//getchar();


		}else
		if (in->dos_free_carrier_stats==mb_look_up_table)
		{
			E=E_free_start;
			sum=0.0;

			gdouble w0=0.0;

			for (e=0;e<E_free_points;e++)
			{

				if (electrons==TRUE)
				{
					rho=(sqrtl(E*Q)/(2.0*PI*PI))*pow((2.0*in->me*m0)/(hbar*hbar),3.0/2.0);
				}else
				{
					rho=(sqrtl(E*Q)/(2.0*PI*PI))*pow((2.0*in->mh*m0)/(hbar*hbar),3.0/2.0);
				}

				f=expl((xpos-E)*Q/(kb*tpos));
				sum+=rho*Q*f*dE_free;

				f=expl((xpos-E-dxr/2.0)*Q/(kb*tpos));
				sum_l+=rho*Q*f*dE_free;

				f=expl((xpos-E+dxr/2.0)*Q/(kb*tpos));
				sum_r+=rho*Q*f*dE_free;
				E+=dE_free;
				//if (electrons==FALSE)
				//{
				//	printf("%Le %Le\n",rho,in->mh);
				//	getchar();
				//}

			}

		}else
		if (in->dos_free_carrier_stats==mb_look_up_table_analytic)
		{
			if (electrons==TRUE)
			{
				sum_l=Nc*exp(((xpos-dxr/2.0)*Q)/(kb*tpos));
				sum=Nc*exp((xpos*Q)/(kb*tpos));
				sum_r=Nc*exp(((xpos+dxr/2.0)*Q)/(kb*tpos));

			}else
			{
				sum_l=Nv*exp(((xpos-dxr/2.0)*Q)/(kb*tpos));
				sum=Nv*exp((xpos*Q)/(kb*tpos));
				sum_r=Nv*exp(((xpos+dxr/2.0)*Q)/(kb*tpos));
			}
		}else
		{
			sum=-1.0;
		}

		//if (electrons==TRUE)
		//{
		//	printf("%Le %Le %Lf\n",in->Nc*exp((xpos*Q)/(kb*tpos)),sum,tpos);
		//}
		//getchar();
		//printf("%Le %Le\n",xpos-dxr/2.0,xpos+dxr/2.0);
		gdouble w0=(3.0/2.0)*sum/((fabs(sum_l-sum_r))/(dxr*Q));

		//printf("%Le %Le\n",w0,(3.0/2.0)*kb*tpos);
		//getchar();
		//if (x==0) w0=(3.0/2.0)*kb*tpos;

		buf[buf_pos++]=sum;
		buf[buf_pos++]=w0;


		for (srh_band=0;srh_band<in->srh_bands;srh_band++)
		{
				buf[buf_pos++]=srh_r1[srh_band];
				buf[buf_pos++]=srh_r2[srh_band];
				buf[buf_pos++]=srh_r3[srh_band];
				buf[buf_pos++]=srh_r4[srh_band];

				buf[buf_pos++]=srh_n[srh_band];
		}

		#ifdef dos_test_stats
			if ((x%100)==0)
			{
				//printf("%d\n",x);
				if (electrons==TRUE)
				{
					freetest=fopen("freetestn.dat","a");
					//fprintf(freetest,"%.20le %.20le %.20le %.20le\n",xpos,sum,in->Nc*exp((xpos*Q)/(kb*tpos)),in->Nv*exp((-(in->Eg+xpos)*Q)/(kb*tpos)));
					//fprintf(freetest,"%Le %Le %Le\n",xpos,in->Nc*exp((xpos*Q)/(kb*tpos)),sum);

					for (srh_band=0;srh_band<in->srh_bands;srh_band++)
					{
						fprintf(freetest,"%Le %Le\n",srh_mid[srh_band],srh_r2[0]);
					}
					fprintf(freetest,"\n");

					fclose(freetest);
				}else
				{
					freetest=fopen("freetestp.dat","a");
					//fprintf(freetest,"%Le %Le %Le\n",xpos,in->Nv*exp((xpos*Q)/(kb*tpos)),sum);

					for (srh_band=0;srh_band<in->srh_bands;srh_band++)
					{
						fprintf(freetest,"%Le %Le\n",srh_mid[srh_band],srh_n[srh_band]);
					}
					fprintf(freetest,"\n");

					fclose(freetest);
				}
			}
		#endif


		#ifdef test_dist
		fprintf(rod,"\n");
		//if (t==0) fprintf(munfile,"%Le %Le\n",n_tot,mu_tot/n_tot);
		#endif


	}
	//getchar();
//getchar();

tpos+=dt;
}

#ifdef test_dist
	fclose(rod);
	fclose(plotbands);
	fclose(plotbandsfree);
	fclose(munfile);
#endif

if (get_dump_status(sim,dump_write_out_band_structure)==TRUE)
{
	if (electrons==TRUE)
	{
		sprintf(name,"%s_dosoutn.dat",confige[mat].dos_name);
	}else
	{
		sprintf(name,"%s_dosoutp.dat",configh[mat].dos_name);
	}

	buffer_dump(sim,name,&dos_out);
	buffer_free(&dos_out);

}

	if (buf_len!=buf_pos)
	{
	ewe(sim,_("Expected dos size is different from generated\n"));
	}

	write_zip_buffer(sim,outfile,buf,buf_len);
	free(buf);

free(xmesh);
free(srh_r1);
free(srh_r2);
free(srh_r3);
free(srh_r4);
free(srh_n);
free(srh_den);
free(srh_dE_sum);
free(srh_read);

free(band_E_mesh);
free(band_i);

free(srh_x);
free(srh_mid);


//if (electrons==TRUE)
//{
//		FILE *test3=fopen("test3.dat","w");
//		for (srh_band=0;srh_band<in->srh_bands;srh_band++)
//		{
//			fprintf(test3,"%d %Le\n",srh_band,srh_den[srh_band]);
//		}
//		fclose(test3);
//		getchar();
//}

return;
}


void gen_dos_fd_gaus_n(struct simulation *sim,int mat)
{

char temp[100];
char path[PATH_MAX];


printf_log(sim,"Electrons.... %s\n",confige[mat].dos_name);

join_path(2, path,get_input_path(sim),"cache");

gpvdm_mkdir(path);

sprintf(temp,"%s_dosn.dat",confige[mat].dos_name);
join_path(3, path,get_input_path(sim),"cache",temp);

gen_do(sim,&confige[mat],&configh[mat],path,TRUE,mat);
}

void gen_dos_fd_gaus_p(struct simulation *sim,int mat)
{
char temp[100];
char path[PATH_MAX];
FILE *out;

printf_log(sim,"Holes.... %s\n",configh[mat].dos_name);

join_path(2, path,get_input_path(sim),"cache");
gpvdm_mkdir(path);

sprintf(temp,"%s_dosp.dat",configh[mat].dos_name);
join_path(3, path,get_input_path(sim),"cache",temp);
gen_do(sim,&configh[mat],&confige[mat],path,FALSE,mat);

join_path(3, path,get_input_path(sim),"cache","mat.inp");

out=fopen(path,"w");
fprintf(out,"#gpvdm_file_type\n");
fprintf(out,"cache\n");
fprintf(out,"#end\n");
fclose(out);
}



void gen_load_dos(struct simulation *sim,int mat,struct epitaxy *my_epitaxy)
{

char file_name[100];
char temp[100];
char full_name[100];
strcpy(confige[mat].dos_name,(my_epitaxy->layer[mat].dos_file));
strcpy(confige[mat].analytical_dos_file_name,(my_epitaxy->lumo_file[mat]));

strcpy(configh[mat].dos_name,(my_epitaxy->layer[mat].dos_file));
strcpy(configh[mat].analytical_dos_file_name,(my_epitaxy->homo_file[mat]));


sprintf(file_name,"%s.inp",confige[mat].dos_name);
join_path(2, full_name,get_input_path(sim),file_name);

struct inp_file inp;
inp_init(sim,&inp);
inp_load(sim,&inp,full_name);
inp_check(sim,&inp,1.24);

inp_search_string(sim,&inp,temp,"#dostype");
confige[mat].dostype=english_to_bin(sim,temp);
configh[mat].dostype=confige[mat].dostype;

inp_search_string(sim,&inp,temp,"#dos_free_carrier_stats");
confige[mat].dos_free_carrier_stats=english_to_bin(sim,temp);
configh[mat].dos_free_carrier_stats=confige[mat].dos_free_carrier_stats;

inp_search_gdouble(sim,&inp,&(confige[mat].Nt),"#Ntrape");
inp_search_gdouble(sim,&inp,&(configh[mat].Nt),"#Ntraph");

confige[mat].Nt=fabs(confige[mat].Nt);
configh[mat].Nt=fabs(configh[mat].Nt);

if (confige[mat].Nt<1e7) confige[mat].Nt=1e7;
if (configh[mat].Nt<1e7) configh[mat].Nt=1e7;


inp_search_gdouble(sim,&inp,&(confige[mat].Et),"#Etrape");
inp_search_gdouble(sim,&inp,&(configh[mat].Et),"#Etraph");

confige[mat].Et=fabs(confige[mat].Et);
configh[mat].Et=fabs(configh[mat].Et);

if (confige[mat].Et<2e-3) confige[mat].Et=2e-3;
if (configh[mat].Et<2e-3) configh[mat].Et=2e-3;

if (confige[mat].Et>200e-3) confige[mat].Et=200e-3;
if (configh[mat].Et>200e-3) configh[mat].Et=200e-3;

inp_search_gdouble(sim,&inp,&(confige[mat].ion_density),"#ion_density");
configh[mat].ion_density=confige[mat].ion_density;

inp_search_gdouble(sim,&inp,&(confige[mat].ion_mobility),"#ion_mobility");
configh[mat].ion_mobility=confige[mat].ion_mobility;

inp_search_gdouble(sim,&inp,&(confige[mat].mu),"#mueffe");
inp_search_gdouble(sim,&inp,&(configh[mat].mu),"#mueffh");

confige[mat].mu=fabs(confige[mat].mu);
configh[mat].mu=fabs(configh[mat].mu);

inp_search_gdouble(sim,&inp,&(confige[mat].epsilonr),"#epsilonr");
confige[mat].epsilonr=fabs(confige[mat].epsilonr);
hard_limit(sim,"#epsilonr",&(confige[mat].epsilonr));

confige[mat].epsilonr=fabs(confige[mat].epsilonr);
configh[mat].epsilonr=fabs(confige[mat].epsilonr);

inp_search_gdouble(sim,&inp,&(confige[mat].doping_start),"#doping_start");
configh[mat].doping_start=confige[mat].doping_start;

inp_search_gdouble(sim,&inp,&(confige[mat].doping_stop),"#doping_stop");
configh[mat].doping_stop=confige[mat].doping_stop;

inp_search_gdouble(sim,&inp,&(confige[mat].Tstart),"#Tstart");
inp_search_gdouble(sim,&inp,&(confige[mat].Tstop),"#Tstop");
inp_search_int(sim,&inp,&(confige[mat].Tsteps),"#Tpoints");

configh[mat].Tstart=confige[mat].Tstart;
configh[mat].Tstop=confige[mat].Tstop;
configh[mat].Tsteps=confige[mat].Tsteps;

inp_search_gdouble(sim,&inp,&(confige[mat].nstart),"#nstart");
inp_search_gdouble(sim,&inp,&(confige[mat].nstop),"#nstop");
inp_search_int(sim,&inp,&(confige[mat].npoints),"#npoints");

inp_search_gdouble(sim,&inp,&(configh[mat].nstart),"#pstart");
inp_search_gdouble(sim,&inp,&(configh[mat].nstop),"#pstop");
inp_search_int(sim,&inp,&(configh[mat].npoints),"#ppoints");
int bands=0;
inp_search_int(sim,&inp,&(bands),"#srh_bands");
confige[mat].srh_bands=bands;
configh[mat].srh_bands=bands;

inp_search_gdouble(sim,&inp,&(confige[mat].srh_start),"#srh_start");
configh[mat].srh_start=confige[mat].srh_start;

inp_search_gdouble(sim,&inp,&(confige[mat].srh_stop),"#srh_stop");
configh[mat].srh_stop=confige[mat].srh_stop;

inp_search_gdouble(sim,&inp,&(confige[mat].srh_sigman),"#srhsigman_e");
confige[mat].srh_sigman=fabs(confige[mat].srh_sigman);

inp_search_gdouble(sim,&inp,&(confige[mat].srh_sigmap),"#srhsigmap_e");
confige[mat].srh_sigmap=fabs(confige[mat].srh_sigmap);

inp_search_gdouble(sim,&inp,&(confige[mat].srh_vth),"#srhvth_e");
confige[mat].srh_vth=fabs(confige[mat].srh_vth);
if (confige[mat].srh_vth<1e2) confige[mat].srh_vth=1e2;

inp_search_gdouble(sim,&inp,&(configh[mat].srh_sigman),"#srhsigman_h");
configh[mat].srh_sigman=fabs(configh[mat].srh_sigman);

inp_search_gdouble(sim,&inp,&(configh[mat].srh_sigmap),"#srhsigmap_h");
configh[mat].srh_sigmap=fabs(configh[mat].srh_sigmap);

inp_search_gdouble(sim,&inp,&(configh[mat].srh_vth),"#srhvth_h");
configh[mat].srh_vth=fabs(configh[mat].srh_vth);
if (configh[mat].srh_vth<1e2) configh[mat].srh_vth=1e2;

inp_search_gdouble(sim,&inp,&(confige[mat].Nc),"#Nc");

inp_search_gdouble(sim,&inp,&(confige[mat].Nv),"#Nv");


inp_search_gdouble(sim,&inp,&(confige[mat].Xi),"#Xi");

inp_search_gdouble(sim,&inp,&(confige[mat].Eg),"#Eg");
confige[mat].Eg=fabs(confige[mat].Eg);
hard_limit(sim,"#Eg",&(confige[mat].Eg));

//confige[mat].Eg=fabs(confige[mat].Eg);
//if (confige[mat].Eg<1.0) configh[mat].Eg=1.0;
//if (confige[mat].Eg>1.8) configh[mat].Eg=1.8;


inp_search_gdouble(sim,&inp,&(confige[mat].B),"#free_to_free_recombination");
confige[mat].B=fabs(confige[mat].B);
configh[mat].B=fabs(confige[mat].B);

//printf("rod >>%Le\n",configh[mat].B);

//gdouble hello;
//inp_search_gdouble(sim,&inp,&(hello),"#free_to_free_recombination");
//printf(">>%Le\n",hello);
//getchar();

configh[mat].B=fabs(confige[mat].B);

int Esteps=0;
int Estep_div=0;

if (bands>0)
{
	inp_search_int(sim,&inp,&(Esteps),"#Esteps");
	Estep_div=(Esteps/bands)*bands;
	//if (Estep_div!=Esteps)
	//{
	//	printf_log(sim,"Esteps wanted= %d, given= %d \n",Esteps,Estep_div);
	//}
}

confige[mat].Esteps=Estep_div;
configh[mat].Esteps=Estep_div;
int dump;
inp_search_int(sim,&inp,&dump,"#dump_band_structure");
set_dump_status(sim,dump_band_structure, dump);

inp_free(sim,&inp);

configh[mat].Xi=confige[mat].Xi;
configh[mat].Eg=confige[mat].Eg;

confige[mat].Nc=fabs(confige[mat].Nc);
confige[mat].Nv=fabs(confige[mat].Nv);
if (confige[mat].Nc<1e16) confige[mat].Nc=1e16;
if (confige[mat].Nv<1e16) confige[mat].Nv=1e16;

/*printf("bing %Lf\n",confige[mat].me);
printf("%Lf\n",confige[mat].mh);
getchar();*/
//confige[mat].m=pow(confige[mat].Nc/2.0,2.0/3.0)*hp*hp/kb/300.0/m0/2.0/PI;
//configh[mat].m=pow(confige[mat].Nv/2.0,2.0/3.0)*hp*hp/kb/300.0/m0/2.0/PI;

//(sqrt(E)/(4.0*PI*PI))*pow((2.0*in->m*m0)/(hbar*hbar),3.0/2.0)

configh[mat].Nc=confige[mat].Nc;
configh[mat].Nv=confige[mat].Nv;

//The states are defined at 300K
confige[mat].me=(powl(confige[mat].Nc/2.0,2.0/3.0)*hbar*hbar*2.0*PI)/(m0*kb*300.0);
confige[mat].mh=(powl(confige[mat].Nv/2.0,2.0/3.0)*hbar*hbar*2.0*PI)/(m0*kb*300.0);

configh[mat].me=confige[mat].me;
configh[mat].mh=confige[mat].mh;

inp_free(sim,&inp);
}

void gen_dos_fd_gaus_fd(struct simulation *sim)
{
char name[200];
char full_name[1000];
int matnumber=0;

struct epitaxy my_epitaxy;
join_path(2, full_name,get_input_path(sim),"epitaxy.inp");
epitaxy_load(sim,&my_epitaxy,full_name);


matnumber=my_epitaxy.electrical_layers;
int file_bandn=FALSE;
int file_bandp=FALSE;
int file_an_lumo=FALSE;
int file_an_homo=FALSE;
int file_dos=FALSE;
int launch_server=FALSE;

FILE *file;
int mat=0;
int problem_with_dos=FALSE;

for (mat=0;mat<matnumber;mat++)
{
	file_bandn=FALSE;
	file_bandp=FALSE;

	file_an_lumo=FALSE;
	file_an_homo=FALSE;
	file_dos=FALSE;

	pick_init(mat);
	gen_load_dos(sim,mat,&my_epitaxy);

	problem_with_dos=FALSE;


	sprintf(name,"%s.inp",my_epitaxy.layer[mat].dos_file);
	join_path(2, full_name,get_input_path(sim),name);

	if (checksum_check(sim,full_name)==FALSE)
	{
		problem_with_dos=TRUE;
	}

	sprintf(name,"%s_dosn.dat",my_epitaxy.layer[mat].dos_file);
	join_path(3, full_name,get_input_path(sim),"cache",name);
	if (isfile(full_name)!=0)
	{
		problem_with_dos=TRUE;
	}

	sprintf(name,"%s_dosp.dat",my_epitaxy.layer[mat].dos_file);
	join_path(3, full_name,get_input_path(sim),"cache",name);
	if (isfile(full_name)!=0)
	{
		problem_with_dos=TRUE;
	}

	if (problem_with_dos==TRUE)
	{
		file_dos=TRUE;
		file_bandn=TRUE;
		file_bandp=TRUE;
		launch_server=TRUE;
	}


	if (confige[mat].dostype==dos_read)
	{
		sprintf(name,"%s_srhbandn.inp",my_epitaxy.layer[mat].dos_file);
		join_path(2, full_name,get_input_path(sim),name);
		if (checksum_check(sim,full_name)==FALSE)
		{
			file_bandn=TRUE;
			launch_server=TRUE;
		}

		sprintf(name,"%s_srhbandp.inp",my_epitaxy.layer[mat].dos_file);
		join_path(2, full_name,get_input_path(sim),name);
		if (checksum_check(sim,full_name)==FALSE)
		{
			file_bandp=TRUE;
			launch_server=TRUE;
		}
	}

	if (confige[mat].dostype==dos_an)
	{
		sprintf(name,"%s.inp",my_epitaxy.lumo_file[mat]);
		join_path(2, full_name,get_input_path(sim),name);
		if (checksum_check(sim,full_name)==FALSE)
		{
			file_an_lumo=TRUE;
			file_bandn=TRUE;
			launch_server=TRUE;
		}

		sprintf(name,"%s.inp",my_epitaxy.homo_file[mat]);
		join_path(2, full_name,get_input_path(sim),name);
		if (checksum_check(sim,full_name)==FALSE)
		{
			file_an_homo=TRUE;
			file_bandp=TRUE;
			launch_server=TRUE;
		}
	}


	if (launch_server==TRUE)
	{
		sprintf(name,"gendosn_%d",mat);
		if (file_bandn==TRUE) server_add_job(sim,name,get_output_path(sim));

		sprintf(name,"gendosp_%d",mat);
		if (file_bandp==TRUE) server_add_job(sim,name,get_output_path(sim));


		pick_dump();

		sprintf(name,"%s.inp",my_epitaxy.layer[mat].dos_file);
		join_path(2, full_name,get_input_path(sim),name);
		if (file_dos==TRUE) checksum_write(sim,full_name);


		if (confige[mat].dostype==dos_read)
		{
			sprintf(name,"%s_srhbandn.inp",my_epitaxy.layer[mat].dos_file);
			safe_file(name);
			join_path(2, full_name,get_input_path(sim),name);
			if (file_bandn==TRUE) checksum_write(sim,full_name);

			sprintf(name,"%s_srhbandp.inp",my_epitaxy.layer[mat].dos_file);
			safe_file(name);
			join_path(2, full_name,get_input_path(sim),name);
			if (file_bandp==TRUE) checksum_write(sim,full_name);
		}

		if (confige[mat].dostype==dos_an)
		{
			sprintf(name,"%s.inp",my_epitaxy.lumo_file[mat]);
			join_path(2, full_name,get_input_path(sim),name);
			if (file_an_lumo==TRUE) checksum_write(sim,full_name);

			sprintf(name,"%s.inp",my_epitaxy.homo_file[mat]);
			join_path(2, full_name,get_input_path(sim),name);
			if (file_an_homo==TRUE) checksum_write(sim,full_name);
		}

	print_jobs(sim);

	server_run_jobs(sim,&(sim->server));
	printf_log(sim,_("Finished generating DoS....\n"));

	}

}



epitaxy_free(sim,&my_epitaxy);
}


