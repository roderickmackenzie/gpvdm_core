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

/** @file dos.c
	@brief Reads in the DoS files but does not generate them, also deals with interpolation.
*/

#include <stdio.h>
#include <stdlib.h>
#include <dos.h>

	#include <zlib.h>
#include "code_ctrl.h"
#include "server.h"
#include "sim.h"
#include "dump.h"
#include "lang.h"
#include "log.h"
#include "cal_path.h"
#include "util.h"

long double get_dos_filled_n(struct device *in)
{
int x=0;
int y=0;
int z=0;

int band;
long double n_tot=0.0;
long double n_tot_den=0.0;

struct dimensions *dim=&(in->ns.dim);

for (z=0;z<dim->zmeshpoints;z++)
{
	for (x=0;x<dim->xmeshpoints;x++)
	{
		for (y=0;y<dim->ymeshpoints;y++)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				n_tot+=in->nt[z][x][y][band];
				n_tot_den+=in->dosn[in->imat[z][x][y]].srh_den[band];
			}
		}
	}
}

n_tot=(n_tot)/(n_tot_den);
return n_tot;
}

long double get_dos_filled_p(struct device *in)
{
int x=0;
int y=0;
int z=0;

int band;
long double p_tot=0.0;
long double p_tot_den=0.0;
struct dimensions *dim=&(in->ns.dim);

for (z=0;z<dim->zmeshpoints;z++)
{
	for (x=0;x<dim->xmeshpoints;x++)
	{
		for (y=0;y<dim->ymeshpoints;y++)
		{
			for (band=0;band<dim->srh_bands;band++)
			{
				p_tot+=in->pt[z][x][y][band];
				p_tot_den+=in->dosp[in->imat[z][x][y]].srh_den[band];
			}
		}
	}
}
p_tot=(p_tot)/(p_tot_den);
return p_tot;
}

long double dos_srh_get_fermi_n(struct device *in,long double n, long double p,int band,int mat,long double T)
{
long double srh_sigman=in->dosn[mat].config.srh_sigman;
long double srh_sigmap=in->dosn[mat].config.srh_sigmap;
long double Nc=in->dosn[mat].config.Nc;
long double Nv=in->dosn[mat].config.Nv;
long double srh_vth=in->dosn[mat].config.srh_vth;
//long double srh_Nt=in->dosn[mat].srh_den[band];
long double srh_en=srh_vth*srh_sigman*Nc*gexp((Q*in->dosn[mat].srh_E[band])/(T*kb));
long double srh_ep=srh_vth*srh_sigmap*Nv*gexp((Q*(-1.0-in->dosn[mat].srh_E[band]))/(T*kb));

long double f=0.0;
f=(n*srh_vth*srh_sigman+srh_ep)/(n*srh_vth*srh_sigman+p*srh_vth*srh_sigmap+srh_en+srh_ep);
long double level=0.0;
level=in->dosn[mat].srh_E[band]-T*kb*logl((1.0/f)-1.0)/Q;
//printf("rod=%Le %Le %Le\n",logl((1.0/f)-1.0),(1.0/f)-1.0,f);
return level;
}

long double dos_srh_get_fermi_p(struct device *in,long double n, long double p,int band,int mat, long double T)
{
long double srh_sigmap=in->dosp[mat].config.srh_sigmap;
long double srh_sigman=in->dosp[mat].config.srh_sigman;
long double Nc=in->dosp[mat].config.Nc;
long double Nv=in->dosp[mat].config.Nv;
long double srh_vth=in->dosp[mat].config.srh_vth;
//long double srh_Nt=in->dosn[mat].srh_den[band];
long double srh_ep=srh_vth*srh_sigmap*Nv*gexp((Q*in->dosp[mat].srh_E[band])/(T*kb));
long double srh_en=srh_vth*srh_sigman*Nc*gexp((Q*(-1.0-in->dosp[mat].srh_E[band]))/(T*kb));
long double f=0.0;
f=(p*srh_vth*srh_sigmap+srh_en)/(p*srh_vth*srh_sigmap+n*srh_vth*srh_sigman+srh_ep+srh_en);
long double level=0.0;
level=in->dosp[mat].srh_E[band]-T*kb*logl((1.0/f)-1.0)/Q;
return level;
}

void dos_init(struct device *in,int mat)
{
in->dosn[mat].used=FALSE;
in->dosp[mat].used=FALSE;
}

void dos_free_now(struct dos *mydos)
{
int ii;
int iii;
if (mydos->used==TRUE)
{
	free(mydos->x);
	free(mydos->t);
	free(mydos->srh_E);
	free(mydos->srh_den);

	for (ii=0;ii<mydos->tlen;ii++)
	{
		for (iii=0;iii<mydos->xlen;iii++)
		{
			free(mydos->srh_r1[ii][iii]);
			free(mydos->srh_r2[ii][iii]);
			free(mydos->srh_r3[ii][iii]);
			free(mydos->srh_r4[ii][iii]);
			free(mydos->srh_c[ii][iii]);

		}

		free(mydos->c[ii]);
		free(mydos->w[ii]);
		free(mydos->srh_r1[ii]);
		free(mydos->srh_r2[ii]);
		free(mydos->srh_r3[ii]);
		free(mydos->srh_r4[ii]);
		free(mydos->srh_c[ii]);
	}

	free(mydos->c);
	free(mydos->w);

	free(mydos->srh_r1);
	free(mydos->srh_r2);
	free(mydos->srh_r3);
	free(mydos->srh_r4);
	free(mydos->srh_c);

	mydos->xlen=0;
	mydos->tlen=0;
	mydos->srh_bands=0;
}

}

void dos_free(struct device *in,int mat)
{
	dos_free_now(&in->dosn[mat]);
	dos_free_now(&in->dosp[mat]);
}

long double get_Nc_free(struct device *in,int mat)
{
return in->dosn[mat].config.Nc;

}

long double get_Nv_free(struct device *in,int mat)
{
return in->dosp[mat].config.Nv;
}

void load_dos_file(struct simulation *sim,struct device *in,struct dos *mydos,char *file)
{
#ifndef dos_bin
long double srh_r1=0.0;
long double srh_r2=0.0;
long double srh_r3=0.0;
long double srh_r4=0.0;
long double srh_c=0.0;
long double w0;
long double n;
#endif
mydos->used=TRUE;
mydos->used=TRUE;
mydos->srh_E=NULL;
mydos->srh_den=NULL;

if (get_dump_status(sim,dump_print_text)==TRUE) printf_log(sim,"%s %s\n",_("Loading file"),file);


	long double *buf;
	int buf_pos=0;
	int buf_len=read_zip_buffer(sim,file,&buf);


	int t=0;
	int x=0;
	int srh_band=0;
	mydos->xlen=(int)buf[buf_pos++];
	mydos->tlen=(int)buf[buf_pos++];
	mydos->srh_bands=(int)buf[buf_pos++];
	mydos->config.epsilonr=buf[buf_pos++];
	mydos->config.doping_start=buf[buf_pos++];
	mydos->config.doping_stop=buf[buf_pos++];
	mydos->config.mu=buf[buf_pos++];
	mydos->config.ion_density=buf[buf_pos++];
	mydos->config.ion_mobility=buf[buf_pos++];
	mydos->config.srh_vth=buf[buf_pos++];
	mydos->config.srh_sigman=buf[buf_pos++];
	mydos->config.srh_sigmap=buf[buf_pos++];
	mydos->config.Nc=buf[buf_pos++];
	mydos->config.Nv=buf[buf_pos++];

	mydos->config.Nt=buf[buf_pos++];
	mydos->config.Eg=buf[buf_pos++];
	mydos->config.Xi=buf[buf_pos++];
	mydos->config.B=buf[buf_pos++];
	mydos->config.dos_free_carrier_stats=(int)buf[buf_pos++];

	long double xsteps=mydos->xlen;
	long double tsteps=mydos->tlen;
	mydos->x=(long double *)malloc(sizeof(long double)*(int)xsteps);
	mydos->t=(long double *)malloc(sizeof(long double)*(int)tsteps);

	if (mydos->srh_bands!=0)
	{
		mydos->srh_E=(long double *)malloc(sizeof(long double)*(int)(mydos->srh_bands));
		mydos->srh_den=(long double *)malloc(sizeof(long double)*(int)(mydos->srh_bands));
	}

	mydos->srh_r1=(long double ***)malloc(sizeof(long double **)*(int)mydos->tlen);
	mydos->srh_r2=(long double ***)malloc(sizeof(long double **)*(int)mydos->tlen);
	mydos->srh_r3=(long double ***)malloc(sizeof(long double **)*(int)mydos->tlen);
	mydos->srh_r4=(long double ***)malloc(sizeof(long double **)*(int)mydos->tlen);
	mydos->srh_c=(long double ***)malloc(sizeof(long double **)*(int)mydos->tlen);

	mydos->w=(long double **)malloc(sizeof(long double **)*(int)mydos->tlen);
	mydos->c=(long double **)malloc(sizeof(long double *)*(int)mydos->tlen);

	for (t=0;t<tsteps;t++)
	{
		mydos->c[t]=(long double *)malloc(sizeof(long double )*(int)xsteps);
		mydos->w[t]=(long double *)malloc(sizeof(long double )*(int)xsteps);
		mydos->srh_r1[t]=(long double **)malloc(sizeof(long double *)*(int)xsteps);
		mydos->srh_r2[t]=(long double **)malloc(sizeof(long double *)*(int)xsteps);
		mydos->srh_r3[t]=(long double **)malloc(sizeof(long double *)*(int)xsteps);
		mydos->srh_r4[t]=(long double **)malloc(sizeof(long double *)*(int)xsteps);
		mydos->srh_c[t]=(long double **)malloc(sizeof(long double *)*(int)xsteps);

		for (x=0;x<xsteps;x++)
		{
			if (mydos->srh_bands!=0)
			{
				mydos->srh_r1[t][x]=(long double *)malloc(sizeof(long double )*(int)(mydos->srh_bands));
				mydos->srh_r2[t][x]=(long double *)malloc(sizeof(long double )*(int)(mydos->srh_bands));
				mydos->srh_r3[t][x]=(long double *)malloc(sizeof(long double )*(int)(mydos->srh_bands));
				mydos->srh_r4[t][x]=(long double *)malloc(sizeof(long double )*(int)(mydos->srh_bands));
				mydos->srh_c[t][x]=(long double *)malloc(sizeof(long double )*(int)(mydos->srh_bands));
			}else
			{
				mydos->srh_r1[t][x]=NULL;
				mydos->srh_r2[t][x]=NULL;
				mydos->srh_r3[t][x]=NULL;
				mydos->srh_r4[t][x]=NULL;
				mydos->srh_c[t][x]=NULL;
			}
		}

	}

	for (x=0;x<xsteps;x++)
	{
		mydos->x[x]=buf[buf_pos++];
		//fprintf(rt,"%le\n",(mydos->x[x]));
	}

	for (t=0;t<tsteps;t++)
	{
		mydos->t[t]=buf[buf_pos++];
		//fprintf(rt,"%le\n",(mydos->t[t]));
	}


	for (srh_band=0;srh_band<(mydos->srh_bands);srh_band++)
	{
		mydos->srh_E[srh_band]=buf[buf_pos++];
		//fprintf(rt,"%le\n",(mydos->srh_E[srh_band]));
	}

	for (srh_band=0;srh_band<(mydos->srh_bands);srh_band++)
	{
		mydos->srh_den[srh_band]=buf[buf_pos++];
		//fprintf(rt,"%le\n",(mydos->srh_den[srh_band]));
	}

	for (t=0;t<tsteps;t++)
	{
		for (x=0;x<xsteps;x++)
		{

			mydos->c[t][x]=buf[buf_pos++];
			mydos->w[t][x]=buf[buf_pos++];
			//fprintf(rt,"%.20le %.20le ",mydos->c[t][x],mydos->w[t][x]);

			for (srh_band=0;srh_band<mydos->srh_bands;srh_band++)
			{
				mydos->srh_r1[t][x][srh_band]=buf[buf_pos++];
				mydos->srh_r2[t][x][srh_band]=buf[buf_pos++];
				mydos->srh_r3[t][x][srh_band]=buf[buf_pos++];
				mydos->srh_r4[t][x][srh_band]=buf[buf_pos++];
				mydos->srh_c[t][x][srh_band]=buf[buf_pos++];

				//fprintf(rt,"%.20le %.20le %.20le %.20le %.20le",mydos->srh_r1[t][x][srh_band],mydos->srh_r2[t][x][srh_band],mydos->srh_r3[t][x][srh_band],mydos->srh_r4[t][x][srh_band],mydos->srh_c[t][x][srh_band]);

			}
		//fprintf(rt,"\n");
		}

	}
	free(buf);

//fclose(rt);

}

long double get_dos_ion_density(struct device *in,int mat)
{
return in->dosn[mat].config.ion_density;
}

long double get_dos_ion_mobility(struct device *in,int mat)
{
return in->dosn[mat].config.ion_mobility;
}


long double get_dos_doping_start(struct device *in,int mat)
{
return in->dosn[mat].config.doping_start;
}

long double get_dos_doping_stop(struct device *in,int mat)
{
return in->dosn[mat].config.doping_stop;
}

long double get_dos_epsilonr(struct device *in,int mat)
{
return in->dosn[mat].config.epsilonr;
}

long double dos_get_band_energy_n(struct device *in,int band,int mat)
{
return in->dosn[mat].srh_E[band];
}

long double dos_get_band_energy_p(struct device *in,int band,int mat)
{
return in->dosp[mat].srh_E[band];
}

long double get_dos_Eg(struct device *in,int mat)
{
return in->dosn[mat].config.Eg;
}

long double get_dos_Xi(struct device *in,int mat)
{
return in->dosn[mat].config.Xi;
}

long double get_dos_B(struct device *in,int mat)
{
return in->dosn[mat].config.B;
}

void load_dos(struct simulation *sim,struct device *in,char *namen, char *namep,int mat)
{
char path[PATH_MAX];

join_path(3, path,sim->root_simulation_path,"cache",namen);
load_dos_file(sim,in,&(in->dosn[mat]),path);

join_path(3, path,sim->root_simulation_path,"cache",namep);
load_dos_file(sim,in,&(in->dosp[mat]),path);
in->ns.dim.srh_bands=in->dosn[mat].srh_bands;

}

long double get_dos_E_n(struct device *in,int band,int mat)
{
return in->dosn[mat].srh_E[band];
}

long double get_dos_E_p(struct device *in,int band,int mat)
{
return in->dosp[mat].srh_E[band];
}

long double get_n_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int r,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;
int t=0;
int x=0;

if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
{
	errors_add(sim,"Electrons asking for %Le but range %Le %Le\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
}


x=hashget(in->dosn[mat].x,in->dosn[mat].xlen,top);
//if (x<0) x=0;


x0=in->dosn[mat].x[x];
x1=in->dosn[mat].x[x+1];
xr=(top-x0)/(x1-x0);

if (in->dosn[mat].tlen>1)
{
	t=hashget(in->dosn[mat].t,in->dosn[mat].tlen,T);
	if (t<0) t=0;

	t0=in->dosn[mat].t[t];
	t1=in->dosn[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}


switch (r)
{
case 1:
	c00=in->dosn[mat].srh_r1[t][x][trap];
	c01=in->dosn[mat].srh_r1[t][x+1][trap];
	c0=c00+xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r1[t+1][x][trap];
		c11=in->dosn[mat].srh_r1[t+1][x+1][trap];
		c1=c10+xr*(c11-c10);
	}
break;

case 2:
	c00=in->dosn[mat].srh_r2[t][x][trap];
	c01=in->dosn[mat].srh_r2[t][x+1][trap];
	c0=c00+xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r2[t+1][x][trap];
		c11=in->dosn[mat].srh_r2[t+1][x+1][trap];
		c1=c10+xr*(c11-c10);
	}
break;

case 3:
	c00=in->dosn[mat].srh_r3[t][x][trap];
	c01=in->dosn[mat].srh_r3[t][x+1][trap];
	c0=c00+xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r3[t+1][x][trap];
		c11=in->dosn[mat].srh_r3[t+1][x+1][trap];
		c1=c10+xr*(c11-c10);
	}
break;

case 4:
	c00=in->dosn[mat].srh_r4[t][x][trap];
	c01=in->dosn[mat].srh_r4[t][x+1][trap];
	c0=c00+xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r4[t+1][x][trap];
		c11=in->dosn[mat].srh_r4[t+1][x+1][trap];
		c1=c10+xr*(c11-c10);
	}
break;
}

c=c0+tr*(c1-c0);
ret=c;

return ret;
}

long double get_p_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int r,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;
int t=0;
int x=0;

if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
{
	errors_add(sim,"Holes asking for %e but range %e %e\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
}

x=hashget(in->dosp[mat].x,in->dosp[mat].xlen,top);
//if (x<0) x=0;

x0=in->dosp[mat].x[x];
x1=in->dosp[mat].x[x+1];
xr=(top-x0)/(x1-x0);


if (in->dosp[mat].tlen>1)
{
	t=hashget(in->dosp[mat].t,in->dosp[mat].tlen,T);
	if (t<0) t=0;

	t0=in->dosp[mat].t[t];
	t1=in->dosp[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}

switch (r)
{
	case 1:
		c00=in->dosp[mat].srh_r1[t][x][trap];
		c01=in->dosp[mat].srh_r1[t][x+1][trap];
		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].srh_r1[t+1][x][trap];
			c11=in->dosp[mat].srh_r1[t+1][x+1][trap];
			c1=c10+xr*(c11-c10);
		}
	break;

	case 2:
		c00=in->dosp[mat].srh_r2[t][x][trap];
		c01=in->dosp[mat].srh_r2[t][x+1][trap];
		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].srh_r2[t+1][x][trap];
			c11=in->dosp[mat].srh_r2[t+1][x+1][trap];
			c1=c10+xr*(c11-c10);
		}
	break;

	case 3:
		c00=in->dosp[mat].srh_r3[t][x][trap];
		c01=in->dosp[mat].srh_r3[t][x+1][trap];
		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].srh_r3[t+1][x][trap];
			c11=in->dosp[mat].srh_r3[t+1][x+1][trap];
			c1=c10+xr*(c11-c10);
		}
	break;

	case 4:
		c00=in->dosp[mat].srh_r4[t][x][trap];
		c01=in->dosp[mat].srh_r4[t][x+1][trap];
		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].srh_r4[t+1][x][trap];
			c11=in->dosp[mat].srh_r4[t+1][x+1][trap];
			c1=c10+xr*(c11-c10);
		}
	break;
}



c=c0+tr*(c1-c0);
ret=c;

return ret;
}



long double get_dn_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int r,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;


int t=0;
int x;

if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
{
	errors_add(sim,"Electrons asking for %Le but range %Le %Le\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
}


x=hashget(in->dosn[mat].x,in->dosn[mat].xlen,top);

x0=in->dosn[mat].x[x];
x1=in->dosn[mat].x[x+1];
xr=1.0/(x1-x0);

if (in->dosn[mat].tlen>1)
{
	t=hashget(in->dosn[mat].t,in->dosn[mat].tlen,T);
	t0=in->dosn[mat].t[t];
	t1=in->dosn[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}

switch (r)
{

case 1:
	c00=in->dosn[mat].srh_r1[t][x][trap];
	c01=in->dosn[mat].srh_r1[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r1[t+1][x][trap];
		c11=in->dosn[mat].srh_r1[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 2:
	c00=in->dosn[mat].srh_r2[t][x][trap];
	c01=in->dosn[mat].srh_r2[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r2[t+1][x][trap];
		c11=in->dosn[mat].srh_r2[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 3:
	c00=in->dosn[mat].srh_r3[t][x][trap];
	c01=in->dosn[mat].srh_r3[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r3[t+1][x][trap];
		c11=in->dosn[mat].srh_r3[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 4:
	c00=in->dosn[mat].srh_r4[t][x][trap];
	c01=in->dosn[mat].srh_r4[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosn[mat].tlen>1)
	{
		c10=in->dosn[mat].srh_r4[t+1][x][trap];
		c11=in->dosn[mat].srh_r4[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;
}

c=c0+tr*(c1-c0);

ret=c;

return ret;
}


long double get_dp_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int r,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;


int t=0;
int x;

if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
{
	errors_add(sim,"Holes asking for %Le but range %Le %Le\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
}


x=hashget(in->dosp[mat].x,in->dosp[mat].xlen,top);

x0=in->dosp[mat].x[x];
x1=in->dosp[mat].x[x+1];
xr=1.0/(x1-x0);

if (in->dosp[mat].tlen>1)
{
	t=hashget(in->dosp[mat].t,in->dosp[mat].tlen,T);
	t0=in->dosp[mat].t[t];
	t1=in->dosp[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}

switch (r)
{
case 1:
	c00=in->dosp[mat].srh_r1[t][x][trap];
	c01=in->dosp[mat].srh_r1[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosp[mat].tlen>1)
	{
		c10=in->dosp[mat].srh_r1[t+1][x][trap];
		c11=in->dosp[mat].srh_r1[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 2:
	c00=in->dosp[mat].srh_r2[t][x][trap];
	c01=in->dosp[mat].srh_r2[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosp[mat].tlen>1)
	{
		c10=in->dosp[mat].srh_r2[t+1][x][trap];
		c11=in->dosp[mat].srh_r2[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 3:
	c00=in->dosp[mat].srh_r3[t][x][trap];
	c01=in->dosp[mat].srh_r3[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosp[mat].tlen>1)
	{
		c10=in->dosp[mat].srh_r3[t+1][x][trap];
		c11=in->dosp[mat].srh_r3[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;

case 4:
	c00=in->dosp[mat].srh_r4[t][x][trap];
	c01=in->dosp[mat].srh_r4[t][x+1][trap];
	c0=xr*(c01-c00);

	if (in->dosp[mat].tlen>1)
	{
		c10=in->dosp[mat].srh_r4[t+1][x][trap];
		c11=in->dosp[mat].srh_r4[t+1][x+1][trap];
		c1=xr*(c11-c10);
	}
break;
}

c=c0+tr*(c1-c0);

ret=c;

return ret;
}



/////////////////////////////////////////////////////trap

long double get_n_pop_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;
int t=0;
int x=0;

if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
{
	errors_add(sim,"Electrons asking for %Le but range %Le %Le\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
}

x=hashget(in->dosn[mat].x,in->dosn[mat].xlen,top);


x0=in->dosn[mat].x[x];
x1=in->dosn[mat].x[x+1];
xr=(top-x0)/(x1-x0);
if (in->dosn[mat].tlen>1)
{
	t=hashget(in->dosn[mat].t,in->dosn[mat].tlen,T);

	t0=in->dosn[mat].t[t];
	t1=in->dosn[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
tr=0.0;
}


c00=in->dosn[mat].srh_c[t][x][trap];
c01=in->dosn[mat].srh_c[t][x+1][trap];
c0=c00+xr*(c01-c00);

if (in->dosn[mat].tlen>1)
{
	c10=in->dosn[mat].srh_c[t+1][x][trap];
	c11=in->dosn[mat].srh_c[t+1][x+1][trap];
	c1=c10+xr*(c11-c10);
}

c=c0+tr*(c1-c0);
ret=c;

return ret;
}

long double get_p_pop_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;
int t=0;
int x=0;

if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
{
	errors_add(sim,"Holes asking for %Le but range %Le %Le\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
}


x=hashget(in->dosp[mat].x,in->dosp[mat].xlen,top);

x0=in->dosp[mat].x[x];
x1=in->dosp[mat].x[x+1];
xr=(top-x0)/(x1-x0);
if (in->dosp[mat].tlen>1)
{
	t=hashget(in->dosp[mat].t,in->dosp[mat].tlen,T);
	t0=in->dosp[mat].t[t];
	t1=in->dosp[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}


c00=in->dosp[mat].srh_c[t][x][trap];
c01=in->dosp[mat].srh_c[t][x+1][trap];
c0=c00+xr*(c01-c00);

if (in->dosp[mat].tlen>1)
{
	c10=in->dosp[mat].srh_c[t+1][x][trap];
	c11=in->dosp[mat].srh_c[t+1][x+1][trap];
	c1=c10+xr*(c11-c10);
}


c=c0+tr*(c1-c0);
ret=c;

return ret;
}



long double get_dn_pop_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap, int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;


int t=0;
int x;
//errors_add(sim,"boo");

if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
{
	errors_add(sim,"Electrons asking for %Le but range %Le %Le\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
}


x=hashget(in->dosn[mat].x,in->dosn[mat].xlen,top);

x0=in->dosn[mat].x[x];
x1=in->dosn[mat].x[x+1];
xr=1.0/(x1-x0);

if (in->dosn[mat].tlen>1)
{
	t=hashget(in->dosn[mat].t,in->dosn[mat].tlen,T);
	t0=in->dosn[mat].t[t];
	t1=in->dosn[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}


c00=in->dosn[mat].srh_c[t][x][trap];
c01=in->dosn[mat].srh_c[t][x+1][trap];
c0=xr*(c01-c00);

if (in->dosn[mat].tlen>1)
{
	c10=in->dosn[mat].srh_c[t+1][x][trap];
	c11=in->dosn[mat].srh_c[t+1][x+1][trap];
	c1=xr*(c11-c10);
}


c=c0+tr*(c1-c0);

ret=c;

return ret;
}


long double get_dp_pop_srh(struct simulation *sim,struct device *in,long double top,long double T,int trap, int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double t0=0.0;
long double t1=0.0;
long double c=0.0;
long double xr=0.0;
long double tr=0.0;
long double c00=0.0;
long double c01=0.0;
long double c10=0.0;
long double c11=0.0;


int t=0;
int x;

if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
{
	errors_add(sim,"Holes asking for %Le but range %Le %Le\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
	//if (get_dump_status(dump_exit_on_dos_error)==TRUE) server_stop_and_exit();
	//exit(0);
}


x=hashget(in->dosp[mat].x,in->dosp[mat].xlen,top);

x0=in->dosp[mat].x[x];
x1=in->dosp[mat].x[x+1];
xr=1.0/(x1-x0);

if (in->dosp[mat].tlen>1)
{
	t=hashget(in->dosp[mat].t,in->dosp[mat].tlen,T);
	t0=in->dosp[mat].t[t];
	t1=in->dosp[mat].t[t+1];
	tr=(T-t0)/(t1-t0);
}else
{
	tr=0.0;
}


c00=in->dosp[mat].srh_c[t][x][trap];
c01=in->dosp[mat].srh_c[t][x+1][trap];
c0=xr*(c01-c00);

if (in->dosp[mat].tlen>1)
{
	c10=in->dosp[mat].srh_c[t+1][x][trap];
	c11=in->dosp[mat].srh_c[t+1][x+1][trap];
	c1=xr*(c11-c10);
}


c=c0+tr*(c1-c0);

ret=c;

return ret;
}


