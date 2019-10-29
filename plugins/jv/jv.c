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

/** @file jv.c
	@brief Simulate JV curve.
*/


#include <sim.h>
#include <exp.h>
#include "jv.h"
#include <dump.h>
#include <dynamic_store.h>
#include "newton_tricks.h"
#include <inp.h>
#include <dat_file.h>
#include <gui_hooks.h>
#include <pl.h>
#include <log.h>
#include <lang.h>
#include <remesh.h>
#include <plot.h>
#include <cal_path.h>
#include <contacts.h>
#include <contacts_vti_store.h>
#include "measure.h"
#include <light_fun.h>
#include <cache.h>

static int unused __attribute__((unused));

int get_step_n(long double step0,long double step_mul,long double V)
{
int n=0;
long double pos=0;
long double dv=fabs(step0);

while(pos<fabs(V))
{
	pos+=dv;
	dv*=step_mul;
	n++;
}
return n;//roundl(log(1.0+(fabs(V)/fabs(step0))*log(step_mul))/log(step_mul));
}

void sim_jv(struct simulation *sim,struct device *in)
{
long double Vapplied=0.0;
int up=TRUE;
struct jv config;
int ittr=0;
gdouble J;
gdouble Pden;
int first=TRUE;
gdouble Vlast;
gdouble Jlast;
gdouble Pdenlast;
gdouble Vexternal;
gdouble V=0.0;

light_solve_and_update(sim,in,&(in->mylight),0.0);


printf_log(sim,_("Running JV simulation\n"));
struct dat_file buf;
buffer_init(&buf);

struct dynamic_store store;
dump_dynamic_init(sim,&store,in);

struct contacts_vti_store contact_store;
dump_contacts_init(sim,in,&contact_store);




char config_file_name[200];

if (find_config_file(sim,config_file_name,get_input_path(sim),in->simmode,"jv")!=0)
{
	ewe(sim,"%s %s %s\n",_("no jv config file found"),get_input_path(sim),in->simmode);
}

printf_log(sim,"%s\n",config_file_name);

jv_load_config(sim,&config,in,config_file_name);

if (config.jv_Rcontact!=-1.0)
{
	in->Rcontact=gfabs(config.jv_Rcontact);
}

if (config.jv_Rshunt!=-1.0)
{
	in->Rshunt=gfabs(config.jv_Rshunt);
}


gdouble Vstop=config.Vstop;
gdouble Vstep=config.Vstep;

struct istruct ivexternal;
inter_init(sim,&ivexternal);

struct istruct jvexternal;
inter_init(sim,&jvexternal);

struct istruct jvavg;
inter_init(sim,&jvavg);

struct istruct charge;
inter_init(sim,&charge);

struct istruct charge_tot;
inter_init(sim,&charge_tot);

struct istruct klist;
inter_init(sim,&klist);

struct istruct lv;
inter_init(sim,&lv);

struct istruct lj;
inter_init(sim,&lj);

//contact_set_active_contact_voltage(sim,in,Vapplied);

//printf("%d\n",in->cache.enabled);
state_cache_enable(sim,in);


if ((in->zmeshpoints>1) || (in->xmeshpoints>1))
{
	contact_set_wanted_active_contact_voltage(sim,in,config.Vstart);
	//contact_set_active_contact_voltage(sim,in,config.Vstart);
	ntricks_auto_ramp_contacts(sim,in);
}else
{
	if (gfabs(config.Vstart-Vapplied)>0.2)
	{
		ramp_externalv(sim,in,Vapplied,config.Vstart);
	}
}


//sim_externalv(in,in->Vapplied);


remesh_reset(in,Vapplied);
//if (in->remesh==TRUE)
//{
//
//}

gdouble sun_orig=light_get_sun(&(in->mylight));
light_set_sun(&(in->mylight),sun_orig*config.jv_light_efficiency);
light_solve_and_update(sim,in,&(in->mylight),0.0);

newton_push_state(in);

newton_set_min_ittr(in,30);

Vapplied=config.Vstart;
contact_set_active_contact_voltage(sim,in,Vapplied);
V=Vapplied;
newton_sim_simple(sim,in);

newton_pop_state(in);
//newton_set_min_ittr(in,0);

//gdouble k_voc=0.0;
gdouble n_voc=0.0;
gdouble r_voc=0.0;
gdouble nsc=0.0;
gdouble n_trap_voc=0.0;
gdouble p_trap_voc=0.0;
gdouble n_free_voc=0.0;
gdouble p_free_voc=0.0;
gdouble np_voc_tot=0.0;
gdouble r_pmax=0.0;
gdouble n_pmax=0.0;
gdouble mue_pmax=0.0;
gdouble muh_pmax=0.0;
long double cal_step=0;

long double n_steps=0.0;
char send_data[200];

struct newton_save_state *ns=&(in->ns);

n_steps=get_step_n(config.Vstep,config.jv_step_mul,config.Vstart);
n_steps+=get_step_n(config.Vstep,config.jv_step_mul,config.Vstop);

in->stop=FALSE;

up=TRUE;
if (config.Vstop<config.Vstart)
{
	up=FALSE;
}

	do
	{
		Vapplied=V;
		contact_set_active_contact_voltage(sim,in,Vapplied);
		newton_sim_simple(sim,in);

		J=get_equiv_J(sim,in);

		Vexternal=get_equiv_V(sim,in);

		sprintf(send_data,"percent:%Lf",(long double)ittr/n_steps);
		gui_send_data(sim,gui_sub,send_data);

		sprintf(send_data,"text:Voltage %.2Lf V/%.2Lf V",V,config.Vstop);
		gui_send_data(sim,gui_sub,send_data);
		if (ittr>0)
		{

			inter_append(&jvexternal,Vexternal,get_equiv_J(sim,in));
			inter_append(&jvavg,V,get_avg_J(in));
			inter_append(&ivexternal,Vexternal,get_equiv_I(sim,in));

		}

		ittr++;
		inter_append(&charge,Vexternal,get_extracted_np(in));
		inter_append(&charge_tot,Vexternal,get_np_tot(in));

		Pden=gfabs(J*Vexternal);

		plot_now(sim,in,"jv.plot");
		stop_start(sim,in);
		dump_dynamic_add_data(sim,&store,in,Vexternal);
		dump_contacts_add_data(sim,in,&contact_store);

		if (get_dump_status(sim,dump_print_converge)==TRUE)
		{
		printf_log(sim," %s=%Lf (%Lf) %s = %Le mA (%Le A/m^2) %Le\n",_("Voltage"),V,Vexternal,_("Current"),get_I(in)/1e-3,J,ns->last_error);
		}

		if (first==FALSE)
		{

			//check if we have crossed 0V
			if ((Vlast<=0)&&(Vexternal>=0.0))
			{
				in->Jsc=Jlast+(J-Jlast)*(0-Vlast)/(V-Vlast);
				nsc=get_extracted_np(in);
				printf_log(sim,"nsc=%Le\n",nsc);
				printf_log(sim,"Jsc = %Le\n",in->Jsc);
			}

			if ((Jlast<=0)&&(J>=0.0))
			{
				in->Voc=Vlast+(Vexternal-Vlast)*(0-Jlast)/(J-Jlast);
				printf_log(sim,"Voc = %Le\n",in->Voc);
				//k_voc=get_avg_recom(in)/pow(get_extracted_np(in),2.0);
				r_voc=get_avg_recom(in);
				n_voc=get_extracted_np(in);
				np_voc_tot=get_total_np(in);
				n_trap_voc=get_n_trapped_charge(in);
				n_free_voc=get_free_n_charge(in);
				p_trap_voc=get_p_trapped_charge(in);
				p_free_voc=get_free_p_charge(in);



			}

			if ((Pden>Pdenlast)&&(Vexternal>0.0)&&(J<0.0))
			{
				in->Pmax=Pden;
				in->Pmax_voltage=Vexternal;
				r_pmax=get_avg_recom(in);
				n_pmax=get_extracted_np(in);
				mue_pmax=get_avg_mue(in);
				muh_pmax=get_avg_muh(in);
			}

			if (up==TRUE)
			{
				if (Vexternal>Vstop)
				{
					printf_log(sim,"%s %Le>%Le\n",_("Stopping because of Vexternal"),Vexternal,Vstop);
					break;
				}
			}else
			{
				if (Vexternal<Vstop)
				{
					printf_log(sim,"%s %Le>%Le\n",_("Stopping because of Vexternal"),Vexternal,Vstop);
					break;
				}
			}


		}




		Jlast=J;
		Vlast=Vexternal;
		Pdenlast=Pden;
		first=FALSE;

		dump_write_to_disk(sim,in);

		long double optical_power_m2=calculate_photon_power_m2(sim,in);
		inter_append(&lv,Vexternal,optical_power_m2);
		//printf("%Le %Le\n",pl_get_light_energy(),in->my_image.extract_eff[lam]);
		inter_append(&lj,J,optical_power_m2);
		//printf("%Le %le\n",get_avg_recom(in),in->my_image.avg_extract_eff);
		V+=Vstep;
		if (config.jv_step_mul>1.0)
		{
			cal_step=get_step_n(config.Vstep,config.jv_step_mul,V);//roundl(log(1.0+(fabs(V)/config.Vstep)*log(config.jv_step_mul))/log(config.jv_step_mul));
			if (cal_step<0)
			{
				cal_step=1.0;
			}
			Vstep=config.Vstep*powl(config.jv_step_mul,cal_step);
			//printf("%Lf %Lf %Lf %Lf\n",Vstep,cal_step,config.Vstep,config.jv_step_mul);
			//getchar();
		}
		//Vstep=fabs(V)(pow(config.jv_step_mul;
		//dialog_set_progress ((in->Vstart+V)/(in->Vstop-in->Vstart));
		if ((up==TRUE)&&(V>Vstop))
		{
			in->stop=TRUE;
		}

		if ((up==FALSE)&&(V<Vstop))
		{
			in->stop=TRUE;
		}

		if (get_equiv_J(sim,in)>config.jv_max_j)
		{
			in->stop=TRUE;
		}

		if (in->stop==TRUE)
		{
			break;
		}

		inter_append(&klist,get_extracted_np(in),get_avg_recom(in)/(pow(get_extracted_np(in),2.0)));

		stop_start(sim,in);


		poll_gui(sim);

		//contacts_detailed_dump(in);
		//contacts_dump(sim,in);
		//dump_contacts_save(sim,in,&contact_store);

		if (config.jv_single_point==TRUE)
		{
			break;
		}
	}while(1);

in->FF=gfabs(in->Pmax/(in->Jsc*in->Voc));

if (get_dump_status(sim,dump_print_text)==TRUE)
{
	printf_log(sim,"Max possible Jsc = %Lf\n",get_max_Jsc(in));
	printf_log(sim,"Voc= %Lf (V)\n",in->Voc);
	printf_log(sim,"Jsc= %Lf (A/m^2)\n",in->Jsc);
	printf_log(sim,"Pmax= %Lf (W/m^2)\n",in->Pmax);
	printf_log(sim,"Pmax %s= %Lf (V)\n",_("Voltage"),in->Pmax_voltage);
	printf_log(sim,"FF= %Lf\n",in->FF*100.0);
	printf_log(sim,"%s= %Lf percent\n",_("Efficiency"),gfabs(in->Pmax/light_get_sun(&(in->mylight))/1000)*100.0);
}

long double added=0.0;
added=get_tot_photons_abs(in);
printf("photon density= %Le\n", added);

if (dumpfiles_should_dump(sim,"sim_info.dat")==0)
{
	FILE *out;
	out=fopena(get_output_path(sim),"sim_info.dat","w");
	fprintf(out,"#ff\n%Lf\n",in->FF);
	fprintf(out,"#pce\n%Lf\n",gfabs(100.0*in->Pmax/(1000.0*light_get_sun(&(in->mylight)))));
	fprintf(out,"#Pmax\n%Lf\n",in->Pmax);
	fprintf(out,"#voc\n%Lf\n",in->Voc);
	fprintf(out,"#voc_tau\n%Le\n",n_voc/r_voc);
	fprintf(out,"#voc_R\n%Le\n",r_voc);
	fprintf(out,"#jv_voc_k\n%Le\n",r_voc/n_voc);
	fprintf(out,"#jv_pmax_n\n%Le\n",n_pmax);
	fprintf(out,"#jv_pmax_tau\n%Le\n",n_pmax/r_pmax);
	fprintf(out,"#jv_pmax_mue\n%Le\n",mue_pmax);
	fprintf(out,"#jv_pmax_muh\n%Le\n",muh_pmax);
	fprintf(out,"#voc_nt\n%Le\n",n_trap_voc);
	fprintf(out,"#voc_pt\n%Le\n",p_trap_voc);
	fprintf(out,"#voc_nf\n%Le\n",n_free_voc);
	fprintf(out,"#voc_pf\n%Le\n",p_free_voc);
	fprintf(out,"#jsc\n%Le\n",in->Jsc);
	fprintf(out,"#jv_jsc_n\n%Le\n",nsc);
	fprintf(out,"#jv_vbi\n%Le\n",in->vbi);
	fprintf(out,"#jv_gen\n%Le\n",get_avg_gen(in));
	fprintf(out,"#voc_np_tot\n%Le\n",np_voc_tot);
	fprintf(out,"#end");
	fclose(out);
}

buffer_malloc(&buf);
buf.y_mul=1.0;
buf.x_mul=1.0;
sprintf(buf.title,"%s - %s",_("Recombination prefactor"),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Recombination prefactor"));
strcpy(buf.y_units,"Volts");
strcpy(buf.data_units,"m^{-6}s^{-1}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_z=0;
buf.x=1;
buf.y=charge.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,klist.x, klist.data, klist.len);
buffer_dump_path(sim,get_output_path(sim),"k.dat",&buf);
buffer_free(&buf);

buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Excess charge density above equilibrium"),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Charge density"));
strcpy(buf.y_units,"Volts");
strcpy(buf.data_units,"m^{-3}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_z=0;
buf.x=1;
buf.y=charge.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,charge.x, charge.data, charge.len);
buffer_dump_path(sim,get_output_path(sim),"charge.dat",&buf);
buffer_free(&buf);

buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Total charge density"),_("App.ed voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Total charge density"));
strcpy(buf.y_units,"Volts");
strcpy(buf.data_units,"m^{-3}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.logscale_z=0;
buf.x=1;
buf.y=charge.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,charge_tot.x, charge_tot.data, charge.len);
buffer_dump_path(sim,get_output_path(sim),"charge_tot.dat",&buf);
buffer_free(&buf);

buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Current density"),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Current density"));
strcpy(buf.y_units,"Volts");
strcpy(buf.data_units,"A m^{-2}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=jvexternal.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,jvexternal.x, jvexternal.data, jvexternal.len);
buffer_dump_path(sim,get_output_path(sim),"jv.dat",&buf);
buffer_free(&buf);

/*
buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Current density"),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Current density"));
strcpy(buf.y_units,_("Volts"));
strcpy(buf.data_units,"A m^{-2}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=jv.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,jv.x, jv.data, jv.len);
buffer_dump_path(sim,get_output_path(sim),"jv_internal.dat",&buf);
buffer_free(&buf);
*/
/*
buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Current density"),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Current density"));
strcpy(buf.y_units,_("Volts"));
strcpy(buf.data_units,"A m^{-2}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=jvavg.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,jvavg.x, jvavg.data, jvavg.len);
buffer_dump_path(sim,get_output_path(sim),"jv_avg.dat",&buf);
buffer_free(&buf);
*/

inter_mul(&jvexternal,in->area);
buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Current "),_("Applied voltage"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,_("Applied Voltage"));
strcpy(buf.data_label,_("Current"));
strcpy(buf.y_units,_("Volts"));
strcpy(buf.data_units,"A");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=jvexternal.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,jvexternal.x, jvexternal.data, jvexternal.len);
buffer_dump_path(sim,get_output_path(sim),"iv.dat",&buf);
buffer_free(&buf);

buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
buf.data_mul=1;
sprintf(buf.title,"%s - %s",_("Voltage"),_("Light flux"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,("Applied Voltage"));
strcpy(buf.data_label,("Light flux"));
strcpy(buf.y_units,"Volts");
strcpy(buf.data_units,"W m^{-2}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=lv.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,lv.x, lv.data, lv.len);
buffer_dump_path(sim,get_output_path(sim),"lv.dat",&buf);
buffer_free(&buf);



buffer_malloc(&buf);
buf.y_mul=1.0;
buf.data_mul=1.0;
sprintf(buf.title,"%s - %s",_("Current density"),_("Light flux"));
strcpy(buf.type,"xy");
strcpy(buf.y_label,("Current density"));
strcpy(buf.data_label,_("Light flux"));
strcpy(buf.y_units,"A m^{-2}");
strcpy(buf.data_units,"W m^{-2}");
buf.logscale_x=0;
buf.logscale_y=0;
buf.x=1;
buf.y=lj.len;
buf.z=1;
buffer_add_info(sim,&buf);
buffer_add_xy_data(sim,&buf,lj.x, lj.data, lj.len);
buffer_dump_path(sim,get_output_path(sim),"lj.dat",&buf);
buffer_free(&buf);

inter_free(&jvexternal);
inter_free(&jvavg);
inter_free(&charge);
inter_free(&ivexternal);
inter_free(&lv);
inter_free(&lj);
inter_free(&klist);

dump_dynamic_save(sim,in,get_output_path(sim),&store);
dump_dynamic_free(sim,in,&store);


dump_contacts_save(sim,in,&contact_store);
dump_contacts_free(sim,in,&contact_store);

light_set_sun(&(in->mylight),sun_orig);
}




void jv_load_config(struct simulation *sim,struct jv* in,struct device *dev, char* config_file_name)
{
	struct inp_file inp;
	inp_init(sim,&inp);
	inp_load_from_path(sim,&inp,get_input_path(sim),config_file_name);
	inp_check(sim,&inp,1.22);
	inp_search_gdouble(sim,&inp,&(in->Vstart),"#Vstart");
	inp_search_gdouble(sim,&inp,&(in->Vstop),"#Vstop");
	inp_search_gdouble(sim,&inp,&(in->Vstep),"#Vstep");
	in->Vstep=fabs(in->Vstep);

	if (in->Vstop<in->Vstart)
	{
		in->Vstep*=-1.0;
	}

	inp_search_gdouble(sim,&inp,&(in->jv_step_mul),"#jv_step_mul");
	inp_search_gdouble(sim,&inp,&(in->jv_light_efficiency),"#jv_light_efficiency");
	inp_search_gdouble(sim,&inp,&(in->jv_max_j),"#jv_max_j");


	inp_search_gdouble(sim,&inp,&(in->jv_Rshunt),"#jv_Rshunt");
	inp_search_gdouble(sim,&inp,&(in->jv_Rcontact),"#jv_Rcontact");
	inp_search_gdouble(sim,&inp,&(in->jv_Rcontact),"#jv_Rcontact");
	in->jv_single_point=inp_search_english(sim,&inp,"#jv_single_point");
	in->jv_light_efficiency=gfabs(in->jv_light_efficiency);
	inp_free(sim,&inp);

}
