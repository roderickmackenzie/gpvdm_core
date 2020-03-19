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
#include <memory.h>
#include <circuit.h>
#include <epitaxy_struct.h>
#include <epitaxy.h>
#include <device_fun.h>
#include <heat.h>
#include <heat_fun.h>


int run_simulation(struct simulation *sim)
{
	struct device cell;
	int enable_electrical=TRUE;

	log_clear(sim);
	struct stat st = {0};

	printf_log(sim,"%s\n",_("Runing simulation"));

	device_init(sim,&cell);
	cache_init(sim);
	cell.onlypos=FALSE;

	dump_init(sim,&cell);

	set_dump_status(sim,dump_stop_plot, FALSE);
	set_dump_status(sim,dump_print_text, TRUE);
	dump_load_config(sim,&cell);

	char temp[PATH_MAX];
	char device_file_path[PATH_MAX];

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

	device_load_math_config(sim,&cell);

	///////////////Mesh and epitaxy//////////////////////////

	join_path(2,device_file_path,get_input_path(sim),"epitaxy.inp");

	epitaxy_load(sim,&(cell.my_epitaxy),device_file_path);

	mesh_obj_load(sim,&(cell.mesh_data));

	dim->zlen=cell.mesh_data.meshdata_z.tot_points;
	dim->xlen=cell.mesh_data.meshdata_x.tot_points;
	dim->ylen=cell.mesh_data.meshdata_y.tot_points;

	mesh_build(sim,&cell);
	device_get_memory(sim,&cell);

	//printf(">>>>>>%Le\n",cell.ns.dim.ymesh[1]);
	//getchar();


	mesh_numerate_points(sim,&cell);

	///////////////////////////////

	load_config(sim,&cell);

	contacts_load(sim,&cell);

	epitaxy_mask(sim,&cell);

	state_cache_init(sim,&cell);
	state_cache_enable(sim,&cell);
	test_complex_solver_init(sim,cell.complex_solver_name);


	device_build_scene(sim,&(cell));

	heat_load_config(sim,&(cell.thermal), &(cell));

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


	printf_log(sim,"%s: %d\n",_("Loading DoS layers"),cell.my_epitaxy.electrical_layers);
	char tempn[200];
	char tempp[200];
	i=0;

	for (i=0;i<cell.my_epitaxy.electrical_layers;i++)
	{
		dos_init(&cell,i);
		printf_log(sim,"%s %d/%d\n",_("Load DoS"),i,cell.my_epitaxy.electrical_layers);
		sprintf(tempn,"%s_dosn.dat",cell.my_epitaxy.layer[i].dos_file);
		sprintf(tempp,"%s_dosp.dat",cell.my_epitaxy.layer[i].dos_file);
		load_dos(sim,&cell,tempn,tempp,i);
	}

	light_init(&cell.mylight);

	solver_init(sim,cell.solver_name);

	if (enable_electrical==TRUE)
	{

		if ((dim->xlen>1)||(dim->zlen>1))
		{
			strcpy(cell.newton_name,"newton_2d");
		}

		newton_init(sim,cell.newton_name);

		//circuit_load(sim);
		circuit_build_device(sim,&(cell.cir),&cell);

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


		for (z=0;z<dim->zlen;z++)
		{
			for (x=0;x<dim->xlen;x++)
			{
				for (y=0;y<dim->ylen;y++)
				{

					depth=cell.ns.dim.ymesh[y]-cell.layer_start[cell.imat[z][x][y]];
					percent=depth/cell.my_epitaxy.layer[cell.imat_epitaxy[z][x][y]].width;
					cell.Nad[z][x][y]=get_dos_doping_start(&cell,cell.imat[z][x][y])+(get_dos_doping_stop(&cell,cell.imat[z][x][y])-get_dos_doping_start(&cell,cell.imat[z][x][y]))*percent;
				}
			}

		}

		update_material_arrays(sim,&cell);


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

		light_load_config(sim,&cell.mylight,&cell);



		//printf("%d %d\n",get_dump_status(sim,dump_optics_verbose), get_dump_status(sim,dump_optics_summary));
		//getchar();
		if ((get_dump_status(sim,dump_optics_verbose)==TRUE) || (get_dump_status(sim,dump_optics_summary)==TRUE))
		{
			light_setup_dump_dir(sim,&cell.mylight);
		}

		light_load_dlls(sim,&cell.mylight);


		circuit_cal_resistance(sim,&(cell.cir),&cell);
		//update_arrays(&cell);

		contacts_force_to_zero(sim,&cell);


		get_initial(sim,&cell,TRUE);

		remesh_shrink(sim,&cell);

		if (cell.math_enable_pos_solver==TRUE)
		{
			for (z=0;z<dim->zlen;z++)
			{
				for (x=0;x<dim->xlen;x++)
				{
					solve_pos(sim,&cell,z,x);
				}
			}
		}




		time_init(sim,&cell);

		matrix_init(&cell.mx);

		solver_realloc(sim,&cell);



		plot_open(sim);


		cell.go_time=FALSE;

		plot_now(sim,&cell,"plot");
		//set_solver_dump_every_matrix(1);

		find_n0(sim,&cell);

		cell.map_start=cell.Ev[0][0][dim->ylen-1];
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
	//		color_cie_load(sim);
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


	if (enable_electrical==FALSE)
	{
		if (cell.thermal.newton_enable_external_thermal==TRUE)
		{
			//printf("Run thermal solver\n");
			heat_solve(sim, &(cell.thermal),&(cell), 0, 0);
		}
	}

	//Clean up
	cache_dump(sim);
	cache_free(sim);
	epitaxy_free(sim,&cell.my_epitaxy);
	contacts_free(sim,&cell);

	for (i=0;i<cell.my_epitaxy.electrical_layers;i++)
	{
		dos_free(&cell,i);
	}

	device_free(sim,&cell);
	mesh_obj_free(sim,&(cell.mesh_data));

	if (enable_electrical==TRUE)
	{
		//color_cie_load(sim);

		plot_close(sim);

		solver_free_memory(sim,&cell);
		newton_interface_free(sim);


	}

	//Free solver dlls
	test_complex_solver_free(sim);
	solver_free(sim);
	printf_log(sim,"%s %i %s\n", _("Solved"), cell.odes, _("Equations"));

	measure(sim);
	dump_clean_cache_files(sim);

return cell.odes;
}

