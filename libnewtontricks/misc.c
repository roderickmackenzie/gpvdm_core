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

/** @file newton_tricks.c
	@brief A collection of helper functions for the newton solver.
*/


#include <stdio.h>
#include <exp.h>
#include "sim.h"
#include "newton_tricks.h"
#include "gui_hooks.h"
#include <plot.h>
#include <cal_path.h>
#include <thermal.h>
#include <contacts.h>
#include <dump.h>
#include <log.h>

static int unused __attribute__((unused));


static gdouble glob_wanted=0.0;

struct newton_math_state math_save_state;

int solve_cur(struct simulation *sim,struct device *in,int z,int x)
{
	int ret=(*sim->dll_solve_cur)(sim,in,z,x);
	if (is_errors(sim)==0)
	{
		errors_dump(sim);
		ewe(sim,"error");
	}
return ret;
}

void newton_push_state(struct device *in)
{
	math_save_state.min_cur_error=in->min_cur_error;
	math_save_state.max_electrical_itt=in->max_electrical_itt;
	math_save_state.newton_min_itt=in->newton_min_itt;
	math_save_state.electrical_clamp=in->electrical_clamp;
	math_save_state.newton_clever_exit=in->newton_clever_exit;
}

void newton_pop_state(struct device *in)
{
	in->min_cur_error=math_save_state.min_cur_error;
	in->max_electrical_itt=math_save_state.max_electrical_itt;
	in->newton_min_itt=math_save_state.newton_min_itt;
	in->electrical_clamp=math_save_state.electrical_clamp;
	in->newton_clever_exit=math_save_state.newton_clever_exit;
}

void ramp_externalv(struct simulation *sim,struct device *in,gdouble from,gdouble to)
{
gdouble V=from;
gdouble dV=0.12;
if ((to-from)<0.0) dV*= -1.0;
printf_log(sim,"dV=%Le\n",dV);
printf_log(sim,"Ramping: from=%Le to=%Le\n",from,to);

if (fabs(to-from)<=fabs(dV)) return;

do
{
	V+=dV;
	if (get_dump_status(sim,dump_print_text)==TRUE) printf_log(sim,"ramp: %Lf %Lf %d\n",V,to,in->kl_in_newton);


	sim_externalv(sim,in,V);

	//plot_now(sim,in,"jv.plot");
	gui_send_data(sim,gui_sub,"pulse");

	if (fabs(V-to)<fabs(dV))
	{
		break;
	}

}while(1);

if (V!=to)
{
	V=to;
	sim_externalv(sim,in,V);
}

return;
}

void ntricks_auto_ramp_contacts(struct simulation *sim,struct device *in)
{
printf_log(sim,"Multidimentional autoramp\n");
int i;
int changed=TRUE;
char send_data[200];
char temp[200];

gdouble Vapplied=0.0;
in->kl_in_newton=FALSE;
solver_realloc(sim,in);

/*if (state_find_vector(sim,in,NULL)==TRUE)
{
	//printf("bing!\n<");
	solve_all(sim,in);
}*/
//getchar();
newton_push_state(in);

in->min_cur_error=1e-8;
in->max_electrical_itt=100;
in->newton_min_itt=3;
in->electrical_clamp=1.0;
in->newton_clever_exit=FALSE;


while (changed==TRUE)
{

	changed=contacts_itterate_to_desired_voltage(sim,in,temp);

	solve_all(sim,in);
	//save_state(sim,in);
	sprintf(send_data,"text:%s",temp);
	gui_send_data(sim,gui_sub,send_data);

	gui_send_data(sim,gui_sub,"pulse");
	poll_gui(sim);
}

newton_pop_state(in);



printf_log(sim,"Finished with multidimentional autoramp\n");
return;

}



void newton_aux_i(struct simulation *sim,struct device *in,gdouble V,gdouble* i,gdouble* didv,gdouble* didphi,gdouble* didxil,gdouble* didxipl,gdouble* didphir,gdouble* didxir,gdouble* didxipr)
{
//<clean>
gdouble i0=*i;
gdouble didv0=*didv;
gdouble didphi0=*didphi;
gdouble didxil0=*didxil;
gdouble didxipl0=*didxipl;
gdouble didphir0=*didphir;
gdouble didxir0=*didxir;
gdouble didxipr0=*didxipr;

*i=glob_wanted-i0-V/in->Rshunt;
*didv= -didv0-1.0/in->Rshunt;
*didphi= -didphi0;
*didxil= -didxil0;
*didxipl= -didxipl0;
*didphir= -didphir0;
*didxir= -didxir0;
*didxipr= -didxipr0;
//</clean>
return;
}


gdouble sim_i(struct simulation *sim,struct device *in,gdouble wantedi)
{
//<clean>
in->kl_in_newton=TRUE;
solver_realloc(sim,in);
glob_wanted=wantedi;
in->newton_aux=&newton_aux_i;
solve_all(sim,in);
//</clean>
return 0.0;
}




void solve_all(struct simulation *sim,struct device *in)
{
int z=0;
int x=0;
int ittr=0;
int cont=TRUE;
struct dimensions *dim=&in->ns.dim;

if (state_search_and_load(sim,in)==TRUE)
{
	in->newton_only_fill_matrix=TRUE;
	solve_cur(sim,in,z,x);
	in->newton_only_fill_matrix=FALSE;
	return;
}


for (z=0;z<dim->zmeshpoints;z++)
{
//	for (x=0;x<in->xmeshpoints;x++)
//	{

		if (in->newton_enable_external_thermal==FALSE)
		{

			solve_cur(sim,in,z,x);
		}else
		{
			do
			{
				solve_cur(sim,in,z,x);
		
				//plot_now(sim,"thermal.plot");
				//getchar();
				solve_thermal(sim,in);
				//plot_now(sim,"thermal.plot");
				//getchar();

				//plot_now(in);
				///getchar();
				if (((in->thermal_conv==TRUE)&&(in->dd_conv==TRUE))||(ittr>10)) cont=FALSE;
				//getchar();
				ittr++;
			}while(cont==TRUE);
		}
//	}
}

save_state(sim,in);

}

void newton_sim_simple(struct simulation  *sim,struct device *in)
{
in->kl_in_newton=FALSE;
solver_realloc(sim,in);

solve_all(sim,in);
}

