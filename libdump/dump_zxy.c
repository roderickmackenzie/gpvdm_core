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

/** @file dump_zxy_slice.c
@brief dump 1 1d slice across the device.
*/

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <sim.h>
#include <dump.h>
#include <dat_file.h>
#include <util.h>
#include <lang.h>
#include <i.h>
#include <exp.h>
#include <dos.h>
#include <memory.h>



void dump_1d_slice(struct simulation *sim,struct device *in,char *out_dir)
{
int x;
int y;
int z;
long double mul_x;
long double mul_y;
long double mul_z;

long double ***temp_3d;
long double **temp_top;
long double **temp_btm;
int band;
char name[100];
char temp[200];
long double Vexternal=get_equiv_V(sim,in);
struct dat_file buf;
buffer_init(&buf);


	struct newton_state *ns=&(in->ns);
	struct dimensions *dim=&in->ns.dim;
	malloc_zxy_gdouble(dim, &temp_3d);

	buffer_add_dir(sim,out_dir);

	cal_J_drift_diffusion(in);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jn_drift.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Drift current density"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron current density (drift)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jn_drift);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jn_diffusion.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Diffusion current density"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron current density (diffusion)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jn_diffusion);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jp_drift.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Drift current density"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Hole current density (drift)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jp_drift);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jp_diffusion.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Diffusion current density"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Hole current density (diffusion)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jp_diffusion);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Ec.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("LUMO"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Ec);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Ev.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("HOMO"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Ev);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Tl.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Lattice temperature"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Temperature"));
	strcpy(buf.data_units,"K");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Tl);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Te.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Electron temperature"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Temperature"));
	strcpy(buf.data_units,"K");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Te);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Th.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Hole temperature"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Temperature"));
	strcpy(buf.data_units,"K");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Te);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Eg.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Band gap"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Eg);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Fn.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Electron quasi Fermi-level"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Fn);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Fp.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Hole quasi Fermi-level"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Fp);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","phi.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Position"),_("Potential"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Potential"));
	strcpy(buf.data_units,"V");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buf.x=dim->xlen;
	buf.y=dim->ylen+2;
	buf.z=dim->zlen;
	buffer_add_info(sim,&buf);

	malloc_zx_gdouble(dim, &temp_top);
	malloc_zx_gdouble(dim, &temp_btm);

	mem_set_zx_gdouble_from_zx_gdouble(dim, temp_top, in->V_y0);
	mem_add_zx_gdouble_from_zx_gdouble(dim, temp_top, in->Vapplied_y0);

	mem_set_zx_gdouble_from_zx_gdouble(dim, temp_btm, in->V_y1);
	mem_add_zx_gdouble_from_zx_gdouble(dim, temp_btm, in->Vapplied_y1);

	buffer_add_3d_device_data_including_boundaries(sim,&buf,in,  ns->phi,temp_top,temp_btm);

	free_zx_gdouble(dim, &temp_btm);
	free_zx_gdouble(dim, &temp_top);

	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","dphi.dat");
	dim_info_to_buf(&buf,dim);
	strcpy(buf.title,_("Potential"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Potential"));
	strcpy(buf.data_units,"V");
	strcpy(buf.section_one,_("Change in 1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);

	zxy_set_gdouble(dim, temp_3d, 0.0);
	three_d_add_gdouble(dim, temp_3d, ns->phi);
	three_d_add_gdouble(dim, temp_3d, in->phi_save);
	buffer_add_3d_data(sim,&buf,dim, temp_3d);

	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jn.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s %Le",_("Current density"),_("Position"),(in->Jn[0][0][0]+in->Jp[0][0][0]+in->Jn[0][0][dim->ylen-1]+in->Jp[0][0][dim->ylen-1])/2.0);
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buf.x=dim->xlen;
	buf.y=dim->ylen+2;
	buf.z=dim->zlen;
	buffer_add_info(sim,&buf);
	//printf("%Le %Le\n",in->Jn_y0[0][0],in->Jn_y1[0][0]);
	//getchar();
	//buffer_add_3d_data(sim,&buf,dim, in->Jn);
	buffer_add_3d_device_data_including_boundaries(sim,&buf,in,  in->Jn,in->Jn_y0,in->Jn_y1);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","Jp.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s %Le",_("Current density"),_("Position"),(in->Jn[0][0][0]+in->Jp[0][0][0]+in->Jn[0][0][dim->ylen-1]+in->Jp[0][0][dim->ylen-1])/2.0);
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buf.x=dim->xlen;
	buf.y=dim->ylen+2;
	buf.z=dim->zlen;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data_including_boundaries(sim,&buf,in,  in->Jp,in->Jp_y0,in->Jp_y1);
	//buffer_add_3d_data(sim,&buf,dim, in->Jp);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jn_x.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Electron current density x"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Electron current density"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jn_x);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","Jp_x.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Electron current density x"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Hole current density"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Jp_x);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jn_plus_Jp.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Total current density (Jn+Jp)"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Total current density (Jn+Jp)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	//sprintf(temp,"%Le %Le\n",dim->ymesh[0]-dim->ymesh[1]/2,get_J_left(in));

	zxy_set_gdouble(dim, temp_3d, 0.0);
	three_d_add_gdouble(dim, temp_3d, in->Jp);
	three_d_add_gdouble(dim, temp_3d, in->Jn);
	buffer_add_3d_data(sim,&buf,dim, temp_3d);

	//sprintf(temp,"%Le %Le\n",dim->ymesh[dim->ylen-1]-dim->ymesh[1]/2,get_J_right(in));
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Jp_drift_plus_diffusion.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Total current density (Jn+Jp)"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Total current density (Jn+Jp)"));
	strcpy(buf.data_units,"A m^{-2}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Transport"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);

	zxy_set_gdouble(dim, temp_3d, 0.0);
	three_d_add_gdouble(dim, temp_3d, in->Jp_drift);
	three_d_add_gdouble(dim, temp_3d, in->Jp_diffusion);
	buffer_add_3d_data(sim,&buf,dim, temp_3d);

	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Fi.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Equlibrium Fermi-level"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Energy"));
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Fi);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","epsilon_r.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Relative permittivity"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Relative permittivity"));
	strcpy(buf.data_units,"au");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Material parameters"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->epsilonr);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","mu_n.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Electron mobility"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.x_label,_("x-position"));
	strcpy(buf.y_label,_("y-position"));
	strcpy(buf.z_label,_("z-position"));
	strcpy(buf.data_label,_("Electron mobility"));
	strcpy(buf.data_units,"m^{2} V^{-1} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Material parameters"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->mun);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","mu_p.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Hole mobility"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Hole mobility"));
	strcpy(buf.data_units,"m^{2} V^{-1} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Material parameters"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->mup);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","mu_n_ft.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s mu0*nf/nall",_("Average electron mobility"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Mobility"));
	strcpy(buf.data_units,"m^{2} V^{-1} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Material parameters"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);

	x=0;
	y=0;
	z=0;

	for (y=0;y<dim->ylen;y++)
	{
		sprintf(temp,"%Le %Le\n",dim->ymesh[y],in->mun[z][x][y]*in->n[z][x][y]/(in->nt_all[z][x][y]+in->n[z][x][y]));
		buffer_add_string(&buf,temp);
	}
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","mu_p_ft.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s mu0*nf/nall",_("Average electron mobility"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Mobility"));
	strcpy(buf.data_units,"m^{2} V^{-1} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Material parameters"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	for (y=0;y<dim->ylen;y++)
	{
		sprintf(temp,"%Le %Le\n",dim->ymesh[y],in->mup[z][x][y]*in->p[z][x][y]/(in->pt_all[z][x][y]+in->p[z][x][y]));
		buffer_add_string(&buf,temp);
	}
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);



	/////////////////////Recombination
	buffer_malloc(&buf);
	sprintf(name,"%s","G_n.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Free electron generation rate"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Generation rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Gn);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","G_p.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Free hole generation rate"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Generation rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Gp);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","R_srh_nfree.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Free electron loss due to SRH processes (including trapping)"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Recombination rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Rn_srh);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","R_srh_pfree.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Free hole loss due to SRH processes (including trapping)"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Recombination rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim, in->Rp_srh);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","R_nfree_to_pfree.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Free electron -> Free hole recombination rate"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Recombination rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->Rfree);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","R_srh_pfree_to_ntrap.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("SRH recombination free holes -> trapped electrons"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->ntrap_to_p);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","R_srh_nfree_to_ptrap.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("SRH recombination free electrons -> trapped holes"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->ptrap_to_n);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","T_srh_pfree_to_ptrap.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("SRH relaxation free holes -> trapped holes"),_("position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim,  in->prelax);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","T_srh_nfree_to_ntrap.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("SRH relaxation free electrons -> trapped electrons"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Rate"));
	strcpy(buf.data_units,"m^{-3} s^{-1}");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_data(sim,&buf,dim, in->nrelax);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	//////////////////////////////////////Charge

	dump_zxy_charge(sim,in,out_dir);


	/////////////////////////
/*
	buffer_malloc(&buf);
	sprintf(name,"%s","fsrhn.dat");
	buf.y_mul=1e9;
	buf.x_mul=1e9;
	buf.z_mul=1e9;
	sprintf(buf.title,"%s - %s",_("Trap fermi level"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.x_label,_("x-position"));
	strcpy(buf.y_label,_("y-position"));
	strcpy(buf.z_label,_("z-position"));
	strcpy(buf.data_label,_("Electron Fermi level"));
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"nm");
	strcpy(buf.z_units,"nm");
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buf.x=dim->xlen;
	buf.y=dim->ylen;
	buf.z=dim->zlen;
	buffer_add_info(sim,&buf);
	for (y=0;y<dim->ylen;y++)
	{
		sprintf(temp,"%Le ",dim->ymesh[y]);
		buffer_add_string(&buf,temp);
		for (band=0;band<in->srh_bands;band++)
		{
			sprintf(temp,"%Le %Le ",in->Fnt[z][x][y][band],-in->phi[z][x][y]-in->Xi[z][x][y]+dos_srh_get_fermi_n(in,in->n[z][x][y], in->p[z][x][y],band,in->imat[z][x][y],in->Te[z][x][y]));
			buffer_add_string(&buf,temp);
		}
		buffer_add_string(&buf,"\n");

	}
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);


	buffer_malloc(&buf);
	sprintf(name,"%s","fsrhh.dat");
	buf.y_mul=1e9;
	buf.x_mul=1e9;
	buf.z_mul=1e9;
	sprintf(buf.title,"%s - %s",_("Trap fermi level"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.x_label,_("x-position"));
	strcpy(buf.y_label,_("y-position"));
	strcpy(buf.z_label,_("z-position"));
	strcpy(buf.data_label,_("Electron Fermi level"));
	strcpy(buf.x_units,"nm");
	strcpy(buf.y_units,"nm");
	strcpy(buf.z_units,"nm");
	strcpy(buf.data_units,"eV");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Recombination"));
	buf.logscale_x=0;
	buf.logscale_y=0;
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buf.x=dim->xlen;
	buf.y=dim->ylen;
	buf.z=dim->zlen;
	buffer_add_info(sim,&buf);
	for (y=0;y<dim->ylen;y++)
	{
		sprintf(temp,"%Le ",dim->ymesh[y]);
		buffer_add_string(&buf,temp);
		for (band=0;band<in->srh_bands;band++)
		{
			sprintf(temp,"%Le %Le ",in->Fpt[z][x][y][band],-in->phi[z][x][y]-in->Xi[z][x][y]-in->Eg[z][x][y]-dos_srh_get_fermi_p(in,in->n[z][x][y], in->p[z][x][y],band,in->imat[z][x][y],in->Th[z][x][y]));
			buffer_add_string(&buf,temp);
		}
		buffer_add_string(&buf,"\n");

	}
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);
*/

	buffer_malloc(&buf);
	sprintf(name,"%s","imat.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Material number"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Number"));
	strcpy(buf.data_units,"au");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Model"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_3d_device_data_int(sim,&buf,in,in->imat);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","passivate_top.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Passivated"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("True/False"));
	strcpy(buf.data_units,"au");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Model"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_2d_device_data_int(sim,&buf,in,in->passivate_y0);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","passivate_btm.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Passivated"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("True/False"));
	strcpy(buf.data_units,"au");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Model"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	buffer_add_2d_device_data_int(sim,&buf,in,in->passivate_y1);
	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);

	buffer_malloc(&buf);
	sprintf(name,"%s","Efield.dat");
	dim_info_to_buf(&buf,dim);
	sprintf(buf.title,"%s - %s",_("Material number"),_("Position"));
	buffer_set_graph_type(&buf,in);
	strcpy(buf.data_label,_("Number"));
	strcpy(buf.data_units,"au");
	strcpy(buf.section_one,_("1D position space output"));
	strcpy(buf.section_two,_("Band structure"));
	buf.time=in->time;
	buf.Vexternal=Vexternal;
	buffer_add_info(sim,&buf);
	gdouble deriv=0.0;
	for (y=dim->ylen-1;y>1;y--)
	{

		deriv= -(ns->phi[z][x][y]-ns->phi[z][x][y-1])/(dim->ymesh[y]-dim->ymesh[y-1]);
		sprintf(temp,"%Le %Le\n",dim->ymesh[y], deriv);
		buffer_add_string(&buf,temp);
	}

	buffer_dump_path(sim,out_dir,name,&buf);
	buffer_free(&buf);





	free_zxy_gdouble(dim, &temp_3d);
}
