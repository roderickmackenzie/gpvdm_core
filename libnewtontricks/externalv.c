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

/** @file newton_externalv.c
	@brief Run the newton solver for an external voltage.
*/


#include <exp.h>
#include "dump.h"
#include "sim.h"
#include <newton_tricks.h>
#include <contacts.h>

static int glob_use_cap=0;


static gdouble glob_wanted_externalv=0.0;




//    ----=----
//    |       |-----
//    |       |    |
//    |       |   ---
//    |       SC
//    |       |   ---
//    |       |    |
//  Vtot      |-----
//    |       R
//    |	      |
//    ----=----


void newton_externv_aux(struct simulation *sim,struct device *in,gdouble V,gdouble* i,gdouble* didv,gdouble* didphi,gdouble* didxil,gdouble* didxipl,gdouble* didphir,gdouble* didxir,gdouble* didxipr)
{
gdouble C=in->C;
gdouble Vapplied_last=0.0;
if (glob_use_cap==FALSE) C=0.0;
gdouble i0=*i;
gdouble didv0=*didv;
gdouble didphi0=*didphi;
gdouble didxil0=*didxil;
gdouble didxipl0=*didxipl;
gdouble didphir0=*didphir;
gdouble didxir0=*didxir;
gdouble didxipr0=*didxipr;
Vapplied_last=contact_get_active_contact_voltage_last(sim,in);
*i=glob_wanted_externalv-(V+((in->Rcontact+in->Rload)*(i0+V/in->Rshunt+C*(V-Vapplied_last)/in->dt)));
*didv= -(1.0+((in->Rcontact+in->Rload)*(didv0+1.0/in->Rshunt+C*(1.0)/in->dt)));
*didphi= -((in->Rcontact+in->Rload)*didphi0);
*didxil= -((in->Rcontact+in->Rload)*didxil0);
*didxipl= -((in->Rcontact+in->Rload)*didxipl0);
*didphir= -((in->Rcontact+in->Rload)*didphir0);
*didxir= -((in->Rcontact+in->Rload)*didxir0);
*didxipr= -((in->Rcontact+in->Rload)*didxipr0);
return;
}


gdouble newton_externv(struct simulation *sim,struct device *in,gdouble Vtot,int usecap)
{
gdouble Vapplied_last=0.0;
gdouble Vapplied=0.0;


gdouble C=in->C;

in->kl_in_newton=TRUE;
solver_realloc(sim,in);
glob_wanted_externalv=Vtot;
glob_use_cap=usecap;
in->newton_aux=&newton_externv_aux;
solve_all(sim,in);
if (glob_use_cap==FALSE) C=0.0;


Vapplied_last=contact_get_active_contact_voltage_last(sim,in);
Vapplied=contact_get_active_contact_voltage(sim,in);
return get_I(in)+Vapplied/in->Rshunt+C*(Vapplied-Vapplied_last)/in->dt;
}

long double newton_externalv_simple(struct simulation *sim,struct device *in,gdouble V)
{
contact_set_active_contact_voltage(sim,in,V);
in->kl_in_newton=FALSE;
solver_realloc(sim,in);
solve_all(sim,in);
return get_I(in);
}





//////////////////////From misc file//////////////////////////////////////////

void newton_aux_externalv(struct simulation *sim,struct device *in,gdouble V,gdouble* i,gdouble* didv,gdouble* didphi,gdouble* didxil,gdouble* didxipl,gdouble* didphir,gdouble* didxir,gdouble* didxipr)
{
gdouble i0=*i;
gdouble didv0=*didv;
gdouble didphi0=*didphi;
gdouble didxil0=*didxil;
gdouble didxipl0=*didxipl;
gdouble didphir0=*didphir;
gdouble didxir0=*didxir;
gdouble didxipr0=*didxipr;

*i=glob_wanted_externalv-(V+in->Rcontact*(i0+V/in->Rshunt));
*didv= -(1.0+in->Rcontact*(didv0+1.0/in->Rshunt));
*didphi= -(in->Rcontact*(didphi0));
*didxil= -(in->Rcontact*(didxil0));
*didxipl= -(in->Rcontact*(didxipl0));
*didphir= -(in->Rcontact*(didphir0));
*didxir= -(in->Rcontact*(didxir0));
*didxipr= -(in->Rcontact*(didxipr0));
return;
}


gdouble sim_externalv_ittr(struct simulation *sim,struct device *in,gdouble wantedv)
{
	//printf("Enter %Le\n",wantedv);
	gdouble Vapplied=0.0;
	gdouble clamp=0.1;
	gdouble step=0.001;
	gdouble e0;
	gdouble e1;
	gdouble i0;
	gdouble i1;
	gdouble deriv;
	gdouble Rs=in->Rcontact;

	in->kl_in_newton=FALSE;
	solver_realloc(sim,in);

	Vapplied=contact_get_active_contact_voltage(sim,in);
	//printf("ok %Le %Le\n",wantedv,Vapplied);

	solve_all(sim,in);
	i0=get_I(in);

	//printf("done %Le %Le\n",wantedv,Vapplied);


	gdouble itot=i0+Vapplied/in->Rshunt;

	e0=fabs(itot*Rs+Vapplied-wantedv);
	Vapplied+=step;
	contact_set_active_contact_voltage(sim,in,Vapplied);

	solve_all(sim,in);

	i1=get_I(in);
	itot=i1+Vapplied/in->Rshunt;

	e1=fabs(itot*Rs+Vapplied-wantedv);

	deriv=(e1-e0)/step;
	step= -e1/deriv;
	//step=step/(1.0+fabs(step/clamp));
	Vapplied+=step;
	contact_set_active_contact_voltage(sim,in,Vapplied);
	int count=0;
	int max=1000;
	do
	{
		e0=e1;
		solve_all(sim,in);
		itot=i1+Vapplied/in->Rshunt;
		e1=fabs(itot*Rs+Vapplied-wantedv);

		deriv=(e1-e0)/step;
		step= -e1/deriv;
		//gdouble clamp=0.01;
		//if (e1<clamp) clamp=e1/100.0;
		//step=step/(1.0+fabs(step/clamp));
		step=step/(1.0+fabs(step/clamp));
		Vapplied+=step;
		//printf("%Le %Le\n",Vapplied,e1);
		contact_set_active_contact_voltage(sim,in,Vapplied);
		if (count>max) break;
		count++;
	}while(e1>1e-8);

	gdouble ret=get_I(in)+Vapplied/in->Rshunt;
	//getchar();
return ret;
}

gdouble sim_externalv(struct simulation *sim,struct device *in,gdouble wantedv)
{


in->kl_in_newton=TRUE;
solver_realloc(sim,in);

glob_wanted_externalv=wantedv;

in->newton_aux=&newton_aux_externalv;

solve_all(sim,in);


return 0.0;
}
