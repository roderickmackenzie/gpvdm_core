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

/** @file sim.c
@brief init sim structure
*/



#include <sim.h>
#include <string.h>

void sim_init(struct simulation *sim)
{
//plotting
	sim->gnuplot= NULL;
	sim->gnuplot_time= NULL;
	sim->converge= NULL;
	sim->tconverge= NULL;
	sim->log_level= -1;
	sim->last_col=0;
	sim->last_nz=0;
	sim->x=NULL;
	sim->Ap=NULL;
	sim->Ai=NULL;
	sim->Ax=NULL;
	sim->b=NULL;
	sim->Tx=NULL;

	//complex solver
	sim->complex_last_col=0;
	sim->complex_last_nz=0;
	sim->complex_x=NULL;
	sim->complex_xz=NULL;
	sim->complex_Ap=NULL;
	sim->complex_Ai=NULL;
	sim->complex_Ax=NULL;
	sim->complex_Az=NULL;

	sim->dll_solve_cur=NULL;
	sim->dll_solver_realloc=NULL;
	sim->dll_solver_free_memory=NULL;
	sim->dll_solver_handle=NULL;
	strcpy(sim->force_sim_mode,"");

	//fit vars
	sim->last_total_error=-1.0;
	sim->gui=FALSE;
	sim->fitting=FALSE;
	sim->dumpfile=NULL;
	sim->dumpfiles=-1;
	sim->html=FALSE;
	sim->bytes_written=0;
	sim->bytes_read=0;
	sim->files_read=0;
	sim->files_written=0;
	sim->T0=300.0;
	sim->D0=243.75;
	sim->n0=1e20;
	sim->x_matrix_offset=0;

	sim->cache_len=-1;
	sim->cache_max=-1;
	sim->cache=NULL;

	sim->mindbustx=FALSE;

	#ifdef dbus
		sim->connection=NULL;
	#endif


	sim->error_log=NULL;
	sim->error_log_size=0;
	sim->error_log_size_max=0;
	sim->errors_logged=0;
	strcpy(sim->server.dbus_finish_signal,"");

	sim->fd_ext_mem=-1;
	sim->fd_ext_memptr=NULL;
	sim->sem_data_for_slave=NULL;
	sim->sem_data_for_master=NULL;
	strcpy(sim->backing_file,"");
	sim->fd_ext_block_size=-1;
}

