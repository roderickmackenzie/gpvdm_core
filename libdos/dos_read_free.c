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
#include <memory.h>



long double get_n_w(struct device *in,long double top,long double T,int mat)
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

int t;
int x;

	if (in->dosn[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(3.0/2.0)*kb*T;
	}else
	{
		#ifdef dos_warn
		if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
		{
			printf("Free electrons asking for %e but range %e %e\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
			//exit(0);
		}
		#endif

		t=search(in->dosn[mat].t,in->dosn[mat].tlen,T);
		x=search(in->dosn[mat].x,in->dosn[mat].xlen,top);


		if (x<0) x=0;
		if (t<0) t=0;

		x0=in->dosn[mat].x[x];
		x1=in->dosn[mat].x[x+1];
		xr=(top-x0)/(x1-x0);

		if (in->dosn[mat].tlen>1)
		{
			t0=in->dosn[mat].t[t];
			t1=in->dosn[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosn[mat].w[t][x];
		c01=in->dosn[mat].w[t][x+1];

		c0=c00+xr*(c01-c00);

		if (in->dosn[mat].tlen>1)
		{
			c10=in->dosn[mat].w[t+1][x];
			c11=in->dosn[mat].w[t+1][x+1];
			c1=c10+xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);
		ret=c;

		//printf(">>rod%Le %Le\n",(3.0/2.0)*kb*T,ret);
		//getchar();
	}

return ret;
}

long double get_p_w(struct device *in,long double top,long double T,int mat)
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

int t;
int x;

	if (in->dosp[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(3.0/2.0)*kb*T;
	}else
	{
		#ifdef dos_warn
		if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
		{
			printf("Free electrons asking for %e but range %e %e\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
			//exit(0);
		}
		#endif

		t=search(in->dosp[mat].t,in->dosp[mat].tlen,T);
		x=search(in->dosp[mat].x,in->dosp[mat].xlen,top);


		if (x<0) x=0;
		if (t<0) t=0;

		x0=in->dosp[mat].x[x];
		x1=in->dosp[mat].x[x+1];
		xr=(top-x0)/(x1-x0);

		if (in->dosp[mat].tlen>1)
		{
			t0=in->dosp[mat].t[t];
			t1=in->dosp[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosp[mat].w[t][x];
		c01=in->dosp[mat].w[t][x+1];

		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].w[t+1][x];
			c11=in->dosp[mat].w[t+1][x+1];
			c1=c10+xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);
		ret=c;

		//printf(">>p rod%Le %Le\n",(3.0/2.0)*kb*T,ret);
		//getchar();
	}

return ret;
}

long double get_dpdT_den(struct device *in,long double top,long double T,int mat)
{
long double ret=0.0;
long double N=in->dosp[mat].config.Nv;
ret= -((top*Q)/kb)*N*gexp((top*Q)/(kb*T))*gpow(T,-2.0);
return ret;
}

long double get_dndT_den(struct device *in,long double top,long double T,int mat)
{
long double ret=0.0;
long double N=in->dosn[mat].config.Nc;
ret= -((top*Q)/kb)*N*gexp((top*Q)/(kb*T))*gpow(T,-2.0);
return ret;
}

long double get_top_from_n(struct device *in,long double n,long double T,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double cll=0.0;
int xx=0;
long double xr=0.0;

int t=0;
int x=0;


	if (in->dosn[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(kb*T/Q)*log((fabs(n))/in->dosn[mat].config.Nc);
	}else
	{
		if (in->dosn[mat].tlen>1)
		{
			for (t=0;t<in->dosn[mat].tlen-1;t++)
			{
				if (in->dosn[mat].t[t]>T) break;
			}
			t--;

		}else
		{
			t=0;
		}

		xx=0;


		for (x=0;x<in->dosn[mat].xlen-1;x++)
		{
			if (in->dosn[mat].c[t][x]>n) break;
		}

		x--;

		if (x<0) x=0;
		if (xx<0) xx=0;
		if (t<0) t=0;

		x0=in->dosn[mat].c[t][x];
		x1=in->dosn[mat].c[t][x+1];
		xr=(n-x0)/(x1-x0);
		if (xr>1) xr=1;
		c0=in->dosn[mat].x[x];
		c1=in->dosn[mat].x[x+1];

		cll=c0+xr*(c1-c0);

		ret=cll;

	}


return ret;
}


long double get_top_from_p(struct device *in,long double p,long double T,int mat)
{
long double ret=0.0;
long double c0=0.0;
long double c1=0.0;
long double x0=0.0;
long double x1=0.0;
long double cll=0.0;
int xx=0;
long double xr=0.0;

long double c;
int t;
int x;

	if (in->dosp[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(kb*T/Q)*log((fabs(p))/in->dosp[mat].config.Nv);
	}else
	{

		for (t=0;t<in->dosp[mat].tlen-1;t++)
		{
		if (in->dosp[mat].t[t]>T) break;
		}
		t--;

		if (t<0) t=0;

		xx=0;

		for (x=0;x<in->dosp[mat].xlen-1;x++)
		{
			if (in->dosp[mat].c[t][x]>p) break;
		}

		x--;

		if (x<0) x=0;
		if (xx<0) xx=0;

		x0=in->dosp[mat].c[t][x];
		x1=in->dosp[mat].c[t][x+1];
		xr=(p-x0)/(x1-x0);
		c0=in->dosp[mat].x[x];
		c1=in->dosp[mat].x[x+1];
		if (xr>1) xr=1;
		cll=c0+xr*(c1-c0);

		c=cll;

		ret=c;
//		printf(">>>>>> %Le %Le\n",ret ,(kb*T/Q)*log((fabs(p))/in->dosp[mat].config.Nv));
//		getchar();
	}

return ret;
}

long double get_n_den(struct device *in,long double top,long double T,int mat)
{
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
	long double ret=0.0;

	if (in->dosn[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb));
	}else
	{

		#ifdef dos_warn
		if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
		{
			errors_add(sim,"Free electrons asking for %Le but range %Le %Le\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
		}
		#endif

		t=search(in->dosn[mat].t,in->dosn[mat].tlen,T);
		x=search(in->dosn[mat].x,in->dosn[mat].xlen,top);

		if (x<0) x=0;
		if (t<0) t=0;

		x0=in->dosn[mat].x[x];
		x1=in->dosn[mat].x[x+1];
		xr=(top-x0)/(x1-x0);

		if (in->dosn[mat].tlen>1)
		{
			t0=in->dosn[mat].t[t];
			t1=in->dosn[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosn[mat].c[t][x];
		c01=in->dosn[mat].c[t][x+1];
		c0=c00+xr*(c01-c00);

		if (in->dosn[mat].tlen>1)
		{
			c10=in->dosn[mat].c[t+1][x];
			c11=in->dosn[mat].c[t+1][x+1];
			c1=c10+xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);

		ret=c;

			//double N=2.0*pow(((2.0*pi*kb*T*m*m0)/(hp*hp)),1.5);
			//long double test=in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb));
			//printf("test = %Le %Le\n",test,ret);
			//getchar();
		//printf("%Lf %Lf %Lf\n",in->dosn[mat].x[x],in->dosn[mat].x[x+1],top);
		//printf(">>%Le %Le\n",in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb)),ret);
		//printf("%Le %Le\n",in->dosn[mat].config.Nc*gexp((Q*in->dosn[mat].x[x])/(T*kb)),in->dosn[mat].c[0][x]);
		//getchar();

		}
	ret=in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb));
	return ret;
}

long double get_p_den(struct device *in,long double top,long double T, int mat)
{
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
	long double ret=0.0;

	if (in->dosp[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=in->dosp[mat].config.Nv*gexp((Q*top)/(T*kb));
	}else
	{

		#ifdef dos_warn
		if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
		{
			errors_add(sim,"Free holes asking for %Le but range %Le %Le\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
		}
		#endif

		t=search(in->dosp[mat].t,in->dosp[mat].tlen,T);
		x=search(in->dosp[mat].x,in->dosp[mat].xlen,top);

		if (x<0) x=0;
		if (t<0) t=0;

		x0=in->dosp[mat].x[x];
		x1=in->dosp[mat].x[x+1];
		xr=(top-x0)/(x1-x0);

		if (in->dosp[mat].tlen>1)
		{
			t0=in->dosp[mat].t[t];
			t1=in->dosp[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosp[mat].c[t][x];
		c01=in->dosp[mat].c[t][x+1];
		c0=c00+xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].c[t+1][x];
			c11=in->dosp[mat].c[t+1][x+1];
			c1=c10+xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);

		ret=c;

			//double N=2.0*pow(((2.0*pi*kb*T*m*m0)/(hp*hp)),1.5);
			//long double test=in->dosp[mat].config.Nv*gexp((Q*top)/(T*kb));
			//printf("test h = %Le %Le\n",test,ret);
			//getchar();
		}
return ret;
}

long double get_n_mu(struct device *in,int mat)
{
return in->dosn[mat].config.mu;
}

long double get_p_mu(struct device *in,int mat)
{
return in->dosp[mat].config.mu;
}

long double get_dn_den(struct device *in,long double top,long double T, int mat)
{
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

	int t;
	int x;

	long double ret=0.0;

	if (in->dosn[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(Q/(T*kb))*in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb));
	}else
	{
		#ifdef dos_warn
		if ((in->dosn[mat].x[0]>top)||(in->dosn[mat].x[in->dosn[mat].xlen-1]<top))
		{
			errors_add(sim,"Free electrons Asking for %e but range %e %e\n",top,in->dosn[mat].x[0],in->dosn[mat].x[in->dosn[mat].xlen-1]);
		}
		#endif

		t=search(in->dosn[mat].t,in->dosn[mat].tlen,T);
		x=search(in->dosn[mat].x,in->dosn[mat].xlen,top);

		x0=in->dosn[mat].x[x];
		x1=in->dosn[mat].x[x+1];
		xr=1.0/(x1-x0);

		if (in->dosn[mat].tlen>1)
		{
			t0=in->dosn[mat].t[t];
			t1=in->dosn[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosn[mat].c[t][x];
		c01=in->dosn[mat].c[t][x+1];
		c0=xr*(c01-c00);

		if (in->dosn[mat].tlen>1)
		{
			c10=in->dosn[mat].c[t+1][x];
			c11=in->dosn[mat].c[t+1][x+1];
			c1=xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);

		ret=c;
	}

ret=(Q/(T*kb))*in->dosn[mat].config.Nc*gexp((Q*top)/(T*kb));
return ret;
}

long double get_dp_den(struct device *in,long double top,long double T, int mat)
{
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

	int t;
	int x;

	long double ret=0.0;

	if (in->dosp[mat].config.dos_free_carrier_stats==mb_equation)
	{
		ret=(Q/(T*kb))*in->dosp[mat].config.Nv*gexp((Q*top)/(T*kb));
	}else
	{
		#ifdef dos_warn
		if ((in->dosp[mat].x[0]>top)||(in->dosp[mat].x[in->dosp[mat].xlen-1]<top))
		{
			errors_add(sim,"Free electrons Asking for %e but range %e %e\n",top,in->dosp[mat].x[0],in->dosp[mat].x[in->dosp[mat].xlen-1]);
		}
		#endif

		t=search(in->dosp[mat].t,in->dosp[mat].tlen,T);
		x=search(in->dosp[mat].x,in->dosp[mat].xlen,top);

		x0=in->dosp[mat].x[x];
		x1=in->dosp[mat].x[x+1];
		xr=1.0/(x1-x0);

		if (in->dosp[mat].tlen>1)
		{
			t0=in->dosp[mat].t[t];
			t1=in->dosp[mat].t[t+1];
			tr=(T-t0)/(t1-t0);
		}else
		{
			tr=0.0;
		}

		c00=in->dosp[mat].c[t][x];
		c01=in->dosp[mat].c[t][x+1];
		c0=xr*(c01-c00);

		if (in->dosp[mat].tlen>1)
		{
			c10=in->dosp[mat].c[t+1][x];
			c11=in->dosp[mat].c[t+1][x+1];
			c1=xr*(c11-c10);
		}

		c=c0+tr*(c1-c0);

		ret=c;

		//printf(">>%Le %Le\n",ret,(Q/(T*kb))*in->dosp[mat].config.Nv*gexp((Q*top)/(T*kb)));
		//getchar();
	}
return ret;
}



