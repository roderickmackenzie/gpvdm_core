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

/** @file gpvdm_const.h
	@brief Physical constants.
*/

#ifndef h_gpvdm_const
#define h_gpvdm_const

//Physical constants
#define STR_MAX	1024
#define epsilon0 (long double)8.85418782e-12			// m-3 kg-1 s4 A2;
#define epsilon0f (float)8.85418782e-12			// m-3 kg-1 s4 A2;

#define mu0f (float)1.25663706e-6			//

#define hp (long double)6.62606896e-34			//J S Wikipeda
#define PI (long double)3.14159265358979323846

#define PIf (float)3.14159265358979323846

#define hbar (long double)(6.62606896e-34/(2.0*PI))	//Calculated
#define kb (long double)1.3806504e-23			//J K-1 Wiki
#define Q (long double)1.602176487e-19			//C Wikipeda
#define m0 (long double)9.10938215e-31 			//Kg Wikipeda
#define cl  (long double)2.99792458e8			//m/s Wikipieda
#define clf  (float)2.99792458e8			//m/s Wikipieda


//SRH constants
#define srh_1	1
#define srh_2	2
#define srh_3	3
#define srh_4	4
#define interface_schottky	 1

//TRUE/FALSE
#define TRUE 1
#define FALSE 0

#define TOP 0
#define BOTTOM 1
#define RIGHT 2
#define LEFT 3

#define ELECTRON 0
#define HOLE 1

#define FIT_SIMPLEX 0
#define FIT_NEWTON 1
#define FIT_BFGS 2

//tpv light
#define tpv_set_light 0
#define tpv_set_voltage 1
#define tpv_mode_laser	0
#define tpv_mode_sun 1

//sim modes
#define IDEAL_DIODE_IDEAL_LOAD 2
#define LOAD 1
#define OPEN_CIRCUIT 0


//dump control
#define dump_optics 1
#define dump_newton 2
#define dump_plot 3
#define dump_stop_plot 4
#define dump_opt_for_fit 5
#define dump_write_converge 6
#define dump_print_text 7
#define dump_exit_on_dos_error 8
#define dump_energy_slice_switch 9
#define dump_zip_files 10
#define dump_lock 11
#define dump_norm_time_to_one 12
#define dump_band_structure 14
#define dump_print_newtonerror 15
#define dump_print_converge 16
#define dump_first_guess 17
#define dump_1d_slices 18
#define dump_print_pos_error 19
#define dump_optics_verbose 20
#define dump_dynamic 22
#define dump_norm_y_axis 24
#define dump_write_out_band_structure 25
#define dump_optics_summary 26
#define dump_optical_probe 27
#define dump_info_text 28
#define dump_optical_probe_spectrum 29
#define dump_ray_trace_map 31
#define dump_file_access_log 32
#define dump_use_cache 34
#define dump_write_headers 37
#define dump_remove_dos_cache 38

//Atempt to put output in files
#define dump_verbosity_key_results	0	
#define dump_verbosity_everything		1
//dos types
#define dos_exp		0
#define dos_an		1
#define dos_fd		2
#define dos_exp_fd 	3
#define dos_read 	5

//free dos types
#define mb_equation 0
#define mb_look_up_table 1
#define fd_look_up_table 2
#define mb_look_up_table_analytic 3

//contact types
#define contact_ohmic 0
#define contact_schottky 1

//Ray tracer
#define ray_run_never	0
#define ray_run_once    1
#define ray_run_step	2
#define RAY_SIM_EDGE	0
#define RAY_VIEWPOINT	1
#define RAY_OBJECT		2
#define ray_emission_single_point 0
#define ray_emission_electrical_mesh 1


#define measure_voltage		0
#define measure_current		1

#define LAYER_ACTIVE 	0
#define LAYER_CONTACT 	1
#define LAYER_OTHER		2


		#include <linux/limits.h>


#endif
