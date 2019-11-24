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

/** @file sim_run.c
@brief run the simulation
*/

#include "util.h"
#include "sim.h"
#include "dos.h"
#include "dump.h"
#include "complex_solver.h"
#include "newton_tricks.h"
#include "log.h"
#include "inp.h"
#include "solver_interface.h"
#include "newton_interface.h"
#include "mesh.h"
#include "remesh.h"
#include "lang.h"
#include <plot.h>
#include "device.h"
#include <cal_path.h>
#include <string.h>
#include <contacts.h>
#include <cache.h>
#include <sys/stat.h>
#include "measure.h"
#include <color.h>
#include <ray_fun.h>
#include <light_fun.h>


int run_simulation(struct simulation *sim)
{
struct device cell;
int enable_electrical=TRUE;

log_clear(sim);
struct stat st = {0};

printf_log(sim,"%s\n",_("Runing simulation"));

device_init(&cell);
cache_init(sim);
cell.onlypos=FALSE;

dump_init(sim,&cell);

set_dump_status(sim,dump_stop_plot, FALSE);
set_dump_status(sim,dump_print_text, TRUE);
dump_load_config(sim,&cell);

char temp[PATH_MAX];

cell.kl_in_newton=FALSE;
struct dimensions *dim=&(cell.ns.dim);

//if (strcmp(outputpath,"")!=0) strcpy(get_output_path(sim),outputpath);

//if (strcmp(inputpath,"")!=0) strcpy(get_input_path(sim),inputpath);


int i;
int z;
int x;
int y;


join_path(2,temp,get_output_path(sim),"error.dat");
remove_file(sim,temp);

join_path(2,temp,get_output_path(sim),"equilibrium");
remove_dir(sim,temp);

join_path(2,temp,get_output_path(sim),"solver");
remove_dir(sim,temp);

if (get_dump_status(sim,dump_newton)==TRUE)
{
	join_path(2,temp,get_output_path(sim),"solver");

	if (stat(temp, &st) == -1)
	{
		mkdir(temp, 0700);
	}
}

//join_path(2,temp,get_output_path(sim),"snapshots");
//remove_dir(sim,temp);
dump_remove_snapshots(sim);

join_path(2,temp,get_output_path(sim),"optics_output");
remove_dir(sim,temp);

join_path(2,temp,get_output_path(sim),"ray_trace");
remove_dir(sim,temp);

join_path(2,temp,get_output_path(sim),"dynamic");
remove_dir(sim,temp);

join_path(2,temp,get_output_path(sim),"frequency");
remove_dir(sim,temp);
load_config(sim,&cell);


cell.pl_enabled=FALSE;
cell.pl_use_experimental_emission_spectra=FALSE;

for (i=0;i<cell.my_epitaxy.layers;i++)
{
	if (cell.my_epitaxy.layer[i].pl_enabled==TRUE)
	{
		cell.pl_enabled=TRUE;
	}

	if (cell.my_epitaxy.layer[i].pl_use_experimental_emission_spectra==TRUE)
	{
		cell.pl_use_experimental_emission_spectra=TRUE;
	}

}


if (strcmp(sim->force_sim_mode,"")!=0)
{
	strcpy(cell.simmode,sim->force_sim_mode);
}


if (strcmp(cell.simmode,"opticalmodel@optics")==0)
{
	enable_electrical=FALSE;
}

if (strcmp(cell.simmode,"fdtd@fdtd")==0)
{
	enable_electrical=FALSE;
}

if (strcmp(cell.simmode,"trace@trace")==0)
{
	enable_electrical=FALSE;
}


if (enable_electrical==TRUE)
{
	solver_init(sim,cell.solver_name);
	newton_init(sim,cell.newton_name);

	printf_log(sim,"%s: %d\n",_("Loading DoS layers"),cell.my_epitaxy.electrical_layers);
	char tempn[100];
	char tempp[100];
	i=0;

	for (i=0;i<cell.my_epitaxy.electrical_layers;i++)
	{
		dos_init(&cell,i);
		printf_log(sim,"%s %d/%d\n",_("Load DoS"),i,cell.my_epitaxy.electrical_layers);
		sprintf(tempn,"%s_dosn.dat",cell.my_epitaxy.dos_file[i]);
		sprintf(tempp,"%s_dosp.dat",cell.my_epitaxy.dos_file[i]);
		load_dos(sim,&cell,tempn,tempp,i);
	}

	device_alloc_traps(&cell);


	if (get_dump_status(sim,dump_write_converge)==TRUE)
	{
	sim->converge = fopena(get_output_path(sim),"converge.dat","w");
	fclose(sim->converge);

	sim->tconverge=fopena(get_output_path(sim),"tconverge.dat","w");
	fclose(sim->tconverge);
	}


	mesh_cal_layer_widths(&cell);

	long double depth=0.0;
	long double percent=0.0;
	long double value=0.0;


	for (z=0;z<dim->zmeshpoints;z++)
	{
		for (x=0;x<dim->xmeshpoints;x++)
		{
			for (y=0;y<dim->ymeshpoints;y++)
			{

				depth=cell.ns.dim.ymesh[y]-cell.layer_start[cell.imat[z][x][y]];
				percent=depth/cell.my_epitaxy.layer[cell.imat_epitaxy[z][x][y]].width;
				cell.Nad[z][x][y]=get_dos_doping_start(&cell,cell.imat[z][x][y])+(get_dos_doping_stop(&cell,cell.imat[z][x][y])-get_dos_doping_start(&cell,cell.imat[z][x][y]))*percent;
			}
		}

	}

	init_mat_arrays(&cell);




	for (z=0;z<dim->zmeshpoints;z++)
	{
		for (x=0;x<dim->xmeshpoints;x++)
		{
			for (y=0;y<dim->ymeshpoints;y++)
			{
				cell.ns.phi[z][x][y]=0.0;
				cell.n[z][x][y]=0.0;
			}
		}
	}

	contacts_load(sim,&cell);


	cell.C=cell.xlen*cell.zlen*epsilon0*cell.epsilonr[0][0][0]/(cell.ylen+cell.other_layers);
	if (get_dump_status(sim,dump_print_text)==TRUE) printf_log(sim,"C=%Le\n",cell.C);
	//printf("%Le\n",cell.C);
	//getchar();
	cell.A=cell.xlen*cell.zlen;
	cell.Vol=cell.xlen*cell.zlen*cell.ylen;

	///////////////////////light model
	char old_model[100];
	gdouble old_Psun=0.0;
	old_Psun=light_get_sun(&cell.mylight);
	light_init(&cell.mylight);
	light_set_dx(&cell.mylight,cell.ns.dim.ymesh[1]-cell.ns.dim.ymesh[0]);

	light_load_config(sim,&cell.mylight,&cell.my_epitaxy);


	//printf("%d %d\n",get_dump_status(sim,dump_optics_verbose), get_dump_status(sim,dump_optics_summary));
	//getchar();
	if ((get_dump_status(sim,dump_optics_verbose)==TRUE) || (get_dump_status(sim,dump_optics_summary)==TRUE))
	{
		light_setup_dump_dir(sim,&cell.mylight);
	}

	light_load_dlls(sim,&cell.mylight);

	///////////////////////

	//update_arrays(&cell);

	contacts_force_to_zero(sim,&cell);

	state_cache_init(sim,&cell);

	get_initial(sim,&cell,TRUE);


	remesh_shrink(sim,&cell);

	if (cell.math_enable_pos_solver==TRUE)
	{
		for (z=0;z<dim->zmeshpoints;z++)
		{
			for (x=0;x<dim->xmeshpoints;x++)
			{
				solve_pos(sim,&cell,z,x);
			}
		}
	}




	time_init(sim,&cell);

	cell.N=0;
	cell.M=0;

	solver_realloc(sim,&cell);



	plot_open(sim);


	cell.go_time=FALSE;

	plot_now(sim,&cell,"plot");
	//set_solver_dump_every_matrix(1);

	find_n0(sim,&cell);

	cell.map_start=cell.Ev[0][0][dim->ymeshpoints-1];
	cell.map_stop=cell.Ec[0][0][0]+1.0;

	//set_solver_dump_every_matrix(0);


	if (cell.onlypos==TRUE)
	{
		join_path(2,temp,get_output_path(sim),"equilibrium");
		dump_1d_slice(sim,&cell,temp);
		cache_dump(sim);
		cache_free(sim);
		epitaxy_free(sim,&cell.my_epitaxy);
		contacts_free(sim,&cell);
		device_free(sim,&cell);
		mesh_obj_free(sim,&(cell.mesh_data));
		color_cie_load(sim);
		return 0;
	}
}



//Load the dll
if (is_domain(cell.simmode)!=0)
{
	char gussed_full_mode[200];
	if (guess_whole_sim_name(sim,gussed_full_mode,get_input_path(sim),cell.simmode)==0)
	{
		printf_log(sim,"I guess we are using running %s\n",gussed_full_mode);
		strcpy(cell.simmode,gussed_full_mode);
	}else
	{
		ewe(sim,"I could not guess which simulation to run from the mode %s\n",cell.simmode);
	}


}



run_electrical_dll(sim,&cell,strextract_domain(cell.simmode));


cache_dump(sim);
cache_free(sim);
epitaxy_free(sim,&cell.my_epitaxy);
contacts_free(sim,&cell);

if (enable_electrical==TRUE)
{

	device_free(sim,&cell);
	mesh_obj_free(sim,&(cell.mesh_data));
	color_cie_load(sim);

	plot_close(sim);

	for (i=0;i<cell.my_epitaxy.electrical_layers;i++)
	{
		dos_free(&cell,i);
	}
	solver_free_memory(sim,&cell);

	newton_interface_free(sim);
	light_free(sim,&cell.mylight);


}

measure(sim);
dump_clean_cache_files(sim);
return cell.odes;
}

