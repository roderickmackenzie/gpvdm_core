//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
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
//

/** @file sim_struct.h
@brief define the sim structure, the sim structure is used to keep all simulation parameters which are physicaly part of the device. Such as dll locations.
*/


#ifndef sim_struct_h
#define sim_struct_h

#include <stdio.h>
#include <server_struct.h>
#include "cache_struct.h"

#ifdef dbus
	#include <dbus/dbus.h>
#endif

	#include <sys/mman.h>
	#include <sys/stat.h>
	#include <fcntl.h>
	#include <unistd.h>

#include <dirent.h>
#include <i_struct.h>

struct dumpfiles_struct
{
char file_name[100];
char path_name[100];
int write_out;
};

struct simulation
{
	//plotting
	FILE *gnuplot;
	FILE *gnuplot_time;
	FILE *converge;
	FILE *tconverge;
	//dump
	int dump_array[100];
	int dumpfiles;
	struct dumpfiles_struct *dumpfile;

	int log_level;
	//paths
	char plugins_path[PATH_MAX];
	char lang_path[PATH_MAX];
	char input_path[PATH_MAX];
	char root_simulation_path[PATH_MAX];
	char output_path[PATH_MAX];
	char share_path[PATH_MAX];
	char exe_path[PATH_MAX];
	char exe_path_dot_dot[PATH_MAX];
	char materials_path[PATH_MAX];
	char cie_color_path[PATH_MAX];
	char emission_path[PATH_MAX];
	char spectra_path[PATH_MAX];
	char home_path[PATH_MAX];
	char shape_path[PATH_MAX];
	char cache_path[PATH_MAX];
	char gpvdm_local_path[PATH_MAX];
	char tmp_path[PATH_MAX];


	//Matrix solver	-	external dll
	int last_col;
	int last_nz;
	double *x;
	int *Ap;
	int *Ai;
	double *Ax;
	double *b;
	double *Tx;
	int x_matrix_offset;

	//Complex matrix solver - external dll
	int c_last_col;
	int c_last_nz;
	double *c_x;
	double *c_xz;
	int *c_Ap;
	int *c_Ai;
	double *c_Ax;
	double *c_Az;
	double *c_b;
	double *c_bz;
	double *c_Tx;
	double *c_Txz;


	//complex solver internal
	int complex_last_col;
	int complex_last_nz;
	double *complex_x;
	double *complex_xz;
	int *complex_Ap;
	int *complex_Ai;
	double *complex_Ax;
	double *complex_Az;

	//Matrix solver dll	- external
	void (*dll_matrix_init)();
	void (*dll_matrix_solve)();
	void (*dll_matrix_solver_free)();
	void *dll_matrix_handle;

	//Complex matrix solver dll	- internal
	void (*dll_complex_matrix_init)();
	void (*dll_complex_matrix_solve)();
	void (*dll_complex_matrix_solver_free)();
	void *dll_complex_matrix_handle;

	//Solve dlls
	int (*dll_solve_cur)();
	int (*dll_solver_realloc)();
	int (*dll_solver_free_memory)();
	void *dll_solver_handle;
	char force_sim_mode[100];

	//Fitting vars
	double last_total_error;
	int fitting;
	struct server_struct server;

	int gui;
	int html;
	long int bytes_written;
	long int bytes_read;
	long int files_read;
	long int files_written;

	long double T0;
	long double D0;
	long double n0;

	int cache_len;
	int cache_max;
	struct cache_item *cache;

	//gui
	int mindbustx;
	#ifdef dbus
		DBusConnection *connection;
	#endif


	struct math_xy cie_x;
	struct math_xy cie_y;
	struct math_xy cie_z;

	char *error_log;
	int error_log_size;
	int error_log_size_max;
	int errors_logged;

	struct lock lock_data;

		int fd_ext_solver;
		int ext_solver_buf_size;
		double *ext_solver_buf;
		char ext_solver_pipe_name[PATH_MAX];
};

#endif

