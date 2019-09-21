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

/** @file initial.c
@brief setup the initial guess for the solvers, this really is just a really bad guess.
*/


#include <stdlib.h>
#include <dos.h>
#include "sim.h"
#include "dump.h"
#include "log.h"
#include <cal_path.h>
#include <dat_file.h>
#include <lang.h>
#include <string.h>
#include <contacts.h>
#include <cal_path.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <contacts.h>

void init_dump(struct simulation *sim,struct device *in)
{
struct dat_file buf;
char out_dir[400];
char name[400];

if (get_dump_status(sim,dump_first_guess)==TRUE)
{
	struct stat st = {0};

	char out_dir[PATH_MAX];
	join_path(2,out_dir,get_output_path(sim),"equilibrium");

	if (stat(out_dir, &st) == -1)
	{
		mkdir(out_dir, 0700);
	}

	strcpy(out_dir,"equilibrium");

	buffer_init(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","init_Fi.dat");
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	sprintf(buf.title,"%s - %s",_("Equilibrium Fermi-level"),_("position"));
	strcpy(buf.type,"xy");
	strcpy(buf.x_label,_("Position"));
	strcpy(buf.y_label,"Fi");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->xmeshpoints;
	buf.y=in->ymeshpoints;
	buf.z=in->zmeshpoints;
	buf.time=in->time;
	buf.Vexternal=0.0;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data(sim,&buf,in,  in->Fi);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","init_Ec.dat");
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	sprintf(buf.title,"%s - %s",_("LUMO"),_("position"));
	strcpy(buf.type,"xy");
	strcpy(buf.x_label,_("Position"));
	strcpy(buf.y_label,"E_{c}");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->xmeshpoints;
	buf.y=in->ymeshpoints;
	buf.z=in->zmeshpoints;
	buf.time=in->time;
	buf.Vexternal=0.0;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data(sim,&buf,in,  in->Ec);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","init_Ev.dat");
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	sprintf(buf.title,"%s - %s",_("HOMO"),_("position"));
	strcpy(buf.type,"xy");
	strcpy(buf.x_label,_("Position"));
	strcpy(buf.y_label,"E_{v}");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->xmeshpoints;
	buf.y=in->ymeshpoints;
	buf.z=in->zmeshpoints;
	buf.time=in->time;
	buf.Vexternal=0.0;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data(sim,&buf,in,  in->Ev);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","init_n.dat");
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	sprintf(buf.title,"%s - %s",_("Electron density"),_("position"));
	strcpy(buf.type,"xy");
	strcpy(buf.x_label,_("Position"));
	strcpy(buf.y_label,"n");
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"m^{-3}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->xmeshpoints;
	buf.y=in->ymeshpoints;
	buf.z=in->zmeshpoints;
	buf.time=in->time;
	buf.Vexternal=0.0;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data(sim,&buf,in,  in->n);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","init_p.dat");
	buf.y_mul=1.0;
	buf.x_mul=1e9;
	sprintf(buf.title,"%s - %s",_("Hole density"),_("position"));
	strcpy(buf.type,"xy");
	strcpy(buf.x_label,_("Position"));
	strcpy(buf.y_label,_("n"));
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"m^{-3}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.x=in->xmeshpoints;
	buf.y=in->ymeshpoints;
	buf.z=in->zmeshpoints;
	buf.time=in->time;
	buf.Vexternal=0.0;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data(sim,&buf,in,  in->p);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);
}
}


void initial_zero(struct simulation *sim,struct device *in)
{
int i;
int z;
int x;
int y;
int band;
in->electron_affinity_right=0.0;
in->electron_affinity_left=0.0;

for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		in->Fi0_top[z][x]=0.0;
		in->Vl[z][x]=0.0;

		in->l_electrons[z][x]=0.0;
		in->l_holes[z][x]=0.0;

		in->r_electrons[z][x]=0.0;
		in->r_holes[z][x]=0.0;

		in->Vr[z][x]=0.0;
	}
}



for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		for (y=0;y<in->ymeshpoints;y++)
		{
			in->Fi[z][x][y]=0.0;

			in->Fn[z][x][y]=0.0;
			in->Fp[z][x][y]=0.0;

			in->phi[z][x][y]=0.0;

			in->x[z][x][y]=0.0;
			in->xp[z][x][y]=0.0;

			in->Ec[z][x][y]=0.0;
			in->Ev[z][x][y]=0.0;

			in->n[z][x][y]=0.0;
			in->p[z][x][y]=0.0;

			in->mun[z][x][y]=0.0;
			in->mup[z][x][y]=0.0;

			for (band=0;band<in->srh_bands;band++)
			{
				in->Fnt[z][x][y][band]= 0.0;
				in->Fpt[z][x][y][band]= 0.0;
				in->xt[z][x][y][band]=0.0;

				in->nt[z][x][y][band]=0.0;
				in->dnt[z][x][y][band]=0.0;

				in->xpt[z][x][y][band]=0.0;
				in->pt[z][x][y][band]=0.0;
				in->dpt[z][x][y][band]=0.0;
			}

		}
	}
}

in->Vbi=0.0;

}

void get_initial(struct simulation *sim,struct device *in)
{
if (strcmp(in->newton_name,"newton_simple")==0)
{
	initial_zero(sim,in);
	return;
}

int i;
int z;
int x;
int y;

gdouble Ef=0.0;
gdouble phi_ramp=0.0;
gdouble Eg=0.0;
gdouble Xi=0.0;
gdouble charge_left=0.0;
gdouble charge_right=0.0;
gdouble top_l=0.0;
gdouble top_r=0.0;


Ef=0.0;
phi_ramp=0.0;
Eg=in->Eg[0][0][0];
Xi=in->Xi[0][0][0];
charge_left=in->lcharge;
charge_right=in->rcharge;
top_l=0.0;
top_r=0.0;
long double left_ref_to_zero=0.0;
long double right_ref_to_zero=0.0;
gdouble delta_phi=0.0;

if (contacts_get_rcharge_type(sim,in)==ELECTRON)
{
	top_r=get_top_from_n(in,charge_right,in->Te[0][0][in->ymeshpoints-1],in->imat[0][0][in->ymeshpoints-1]);
	in->electron_affinity_right= -in->Xi[0][0][in->ymeshpoints-1]+top_r;
	right_ref_to_zero=top_r-in->Xi[0][0][in->ymeshpoints-1];
}else
{
	top_r= get_top_from_p(in,charge_right,in->Te[0][0][in->ymeshpoints-1],in->imat[0][0][in->ymeshpoints-1]);
	in->electron_affinity_right= -in->Xi[0][0][in->ymeshpoints-1]-in->Eg[0][0][in->ymeshpoints-1]-get_top_from_p(in,charge_right,in->Te[0][0][in->ymeshpoints-1],in->imat[0][0][in->ymeshpoints-1]);
	right_ref_to_zero=-(in->Eg[0][0][in->ymeshpoints-1]+top_r)-in->Xi[0][0][in->ymeshpoints-1];
}

//in->electron_affinity_left= -in->Xi[0][0][0]-in->Eg[0][0][0]-get_top_from_p(in,charge_left,in->Te[0][0][0],in->imat[0][0][0]);
//in->electron_affinity_left= -in->Xi[0][0][0]+top_l;
in->electron_affinity_left=0.0;
int c=0;
int type=0;
for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		c=in->n_contact_l[z][x];

		if (c!=-1)
		{
			if (in->contacts[c].charge_type==HOLE)
			{
				top_l=get_top_from_p(in,in->contacts[c].np,in->Te[z][x][0],in->imat[z][x][0]);
				in->Fi0_top[z][x]= -(top_l+Xi+Eg);
			}else
			{
				top_l= get_top_from_n(in,in->contacts[c].np,in->Te[z][x][0],in->imat[z][x][0]);
				in->Fi0_top[z][x]= -Xi+top_l;
			}
		}else
		{		//No contact
			top_l=get_top_from_p(in,1e15,in->Te[z][x][0],in->imat[z][x][0]);
			in->Fi0_top[z][x]= -(top_l+Xi+Eg);
		}

		in->Vl[z][x]=in->Fi0_top[z][x]-in->Fi0_top[0][0];		//Everything is referenced to the [0][0] point.
		//printf("%Le\n",in->Vl[z][x]);
	}
}
//getchar();
left_ref_to_zero=in->Fi0_top[0][0];
Ef=in->Fi0_top[0][0];


delta_phi=right_ref_to_zero-left_ref_to_zero;


in->vbi=delta_phi;
if (get_dump_status(sim,dump_print_text)==TRUE)
{
	printf_log(sim,"delta=%Le\n",delta_phi);
	printf_log(sim,">>>>top_l= %Le\n",top_l+Eg);
	printf_log(sim,">>>>top_r= %Le\n",-top_r);
	printf_log(sim,"left= %Le right = %Le  %Le %Le\n",in->electron_affinity_left,in->electron_affinity_right,in->electron_affinity_right-in->electron_affinity_left,delta_phi);
	printf_log(sim,"%Le %Le %Le %Le %Le\n",top_l,top_r,Eg,delta_phi,in->phi[0][0][0]);
}




//printf("total %Le\n",(-in->Xi[0][0][0]-in->phi[0][0][0]-Eg)-Ef);

printf(">>rod>>%Le\n",Ef-(-in->Xi[0][0][0]-in->phi[0][0][0]));

gdouble Rp=get_p_den(in,(-in->Xi[0][0][in->ymeshpoints-1]-delta_phi-Eg)-Ef,in->Th[0][0][in->ymeshpoints-1],in->imat[0][0][in->ymeshpoints-1]);
gdouble Rn=get_n_den(in,Ef-(-in->Xi[0][0][in->ymeshpoints-1]-delta_phi),in->Te[0][0][in->ymeshpoints-1],in->imat[0][0][in->ymeshpoints-1]);

for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		in->l_electrons[z][x]=get_n_den(in,Ef-(-in->Xi[z][x][0]-in->Vl[z][x]),in->Te[z][x][0],in->imat[z][x][0]);
		in->l_holes[z][x]=get_p_den(in,(-in->Xi[z][x][0]-in->Vl[z][x]-in->Eg[z][x][0])-Ef,in->Th[z][x][0],in->imat[z][x][0]);;
		printf_log(sim,"Left (%d,%d)  p=%Le n=%Le\n",z,x,in->l_holes[z][x],in->l_electrons[z][x]);

		in->r_electrons[z][x]=Rn;
		in->r_holes[z][x]=Rp;

		in->Vr[z][x]=delta_phi;

	}
}
printf_log(sim,"Rp = %Le\n",Rp);
printf_log(sim,"Rn = %Le\n",Rn);

int band;
for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		for (y=0;y<in->ymeshpoints;y++)
		{
			phi_ramp=delta_phi*(in->ymesh[y]/in->ymesh[in->ymeshpoints-1]);
			//printf("%ld %ld %ld %Le\n",x,y,z,phi_ramp);
			in->Fi[z][x][y]=Ef;

			in->Fn[z][x][y]=Ef;
			in->Fp[z][x][y]=Ef;

			in->phi[z][x][y]=phi_ramp;

			in->x[z][x][y]=in->phi[z][x][y]+in->Fn[z][x][y];
			in->xp[z][x][y]= -(in->phi[z][x][y]+in->Fp[z][x][y]);

			in->Ec[z][x][y]= -in->phi[z][x][y]-in->Xi[z][x][y];
			if (in->Ec[z][x][y]<in->Fi[z][x][y])
			{
				in->phi[z][x][y]= -(in->Fi[z][x][y]+in->Xi[z][x][y]);
				in->Ec[z][x][y]= -in->phi[z][x][y]-in->Xi[z][x][y];
			}

			in->Ev[z][x][y]= -in->phi[z][x][y]-in->Xi[z][x][y]-in->Eg[z][x][y];
			if (in->Ev[z][x][y]>in->Fi[z][x][y])
			{
				in->phi[z][x][y]= -(in->Fi[z][x][y]+in->Xi[z][x][y]+in->Eg[z][x][y]);
				in->Ev[z][x][y]= -in->phi[z][x][y]-in->Xi[z][x][y]-in->Eg[z][x][y];

				in->Ec[z][x][y]= -in->phi[z][x][y]-in->Xi[z][x][y];
			}


			gdouble t=in->Fi[z][x][y]-in->Ec[z][x][y];
			gdouble tp=in->Ev[z][x][y]-in->Fi[z][x][y];

			in->n[z][x][y]=in->Nc[z][x][y]*exp(((t)*Q)/(kb*in->Te[z][x][y]));
			in->p[z][x][y]=in->Nv[z][x][y]*exp(((tp)*Q)/(kb*in->Th[z][x][y]));

			in->mun[z][x][y]=get_n_mu(in,in->imat[z][x][y]);
			in->mup[z][x][y]=get_p_mu(in,in->imat[z][x][y]);

			for (band=0;band<in->srh_bands;band++)
			{
				in->Fnt[z][x][y][band]= Ef;//-in->phi[z][x][y]-in->Xi[z][x][y]+dos_srh_get_fermi_n(in,in->n[z][x][y], in->p[z][x][y],band,in->imat[z][x][y],in->Te[z][x][y]);
				//printf("d %ld %Le\n",band,dos_srh_get_fermi_n(in,in->n[z][x][y], in->p[z][x][y],band,in->imat[z][x][y],in->Te[z][x][y]));
				in->Fpt[z][x][y][band]= Ef;//-in->phi[z][x][y]-in->Xi[z][x][y]-in->Eg[z][x][y]-dos_srh_get_fermi_p(in,in->n[z][x][y], in->p[z][x][y],band,in->imat[z][x][y],in->Th[z][x][y]);
				in->xt[z][x][y][band]=in->phi[z][x][y]+in->Fnt[z][x][y][band];

				in->nt[z][x][y][band]=get_n_pop_srh(sim,in,in->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,in->imat[z][x][y]);
				in->dnt[z][x][y][band]=get_dn_pop_srh(sim,in,in->xt[z][x][y][band]+in->tt[z][x][y],in->Te[z][x][y],band,in->imat[z][x][y]);


				in->xpt[z][x][y][band]= -(in->phi[z][x][y]+in->Fpt[z][x][y][band]);
				in->pt[z][x][y][band]=get_p_pop_srh(sim,in,in->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,in->imat[z][x][y]);
				in->dpt[z][x][y][band]=get_dp_pop_srh(sim,in,in->xpt[z][x][y][band]-in->tpt[z][x][y],in->Th[z][x][y],band,in->imat[z][x][y]);
			}

		}
	}
}

in->Vbi=delta_phi;
init_dump(sim,in);

contacts_passivate(sim,in);
if (in->stoppoint==1) exit(0);
return;
}

