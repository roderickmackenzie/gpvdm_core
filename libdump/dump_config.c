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

/** @file dump_config.c
@brief read the dump config file
*/


#include "sim.h"
#include "dump.h"
#include "inp.h"
#include "log.h"
#include "dump.h"
#include <cal_path.h>

void dump_load_config(struct simulation* sim,struct device *in)
{
	int dump;
	struct inp_file inp;
	inp_init(sim,&inp);

	inp_load_from_path(sim,&inp,get_input_path(sim),"dump.inp");

	inp_check(sim,&inp,1.55);

	dump=inp_search_english(sim,&inp,"#plot");
	set_dump_status(sim,dump_plot,dump);

	dump=inp_search_english(sim,&inp,"#newton_dump");
	set_dump_status(sim,dump_newton, dump);

	dump=inp_search_english(sim,&inp,"#dump_dynamic");
	set_dump_status(sim,dump_dynamic, dump);

	dump=inp_search_english(sim,&inp,"#dump_fx");
	set_dump_status(sim,dump_fx, dump);

	in->stop_start=inp_search_english(sim,&inp,"#startstop");

	in->dumpitdos=inp_search_english(sim,&inp,"#dumpitdos");

	inp_search_string(sim,&inp,in->plot_file,"#plotfile");

	inp_search_gdouble(sim,&inp,&(in->start_stop_time),"#start_stop_time");

	dump=inp_search_english(sim,&inp,"#dump_optics");
	set_dump_status(sim,dump_optics, dump);

	dump=inp_search_english(sim,&inp,"#dump_optics_verbose");
	set_dump_status(sim,dump_optics_verbose, dump);

	dump=inp_search_english(sim,&inp,"#dump_norm_time_to_one");
	set_dump_status(sim,dump_norm_time_to_one, dump);

	dump=inp_search_english(sim,&inp,"#dump_norm_y_axis");
	set_dump_status(sim,dump_norm_y_axis, dump);

	set_dump_status(sim,dump_energy_slice_switch, TRUE);

	inp_search_int(sim,&inp,&(in->dump_energy_slice_xpos),"#dump_energy_slice_xpos");
	if (in->dump_energy_slice_xpos<0)
	{
			set_dump_status(sim,dump_energy_slice_switch, FALSE);
	}
	if (in->dump_energy_slice_xpos>=in->xmeshpoints) in->dump_energy_slice_xpos=0;

	inp_search_int(sim,&inp,&(in->dump_energy_slice_ypos),"#dump_energy_slice_ypos");
	if (in->dump_energy_slice_xpos<0)
	{
			set_dump_status(sim,dump_energy_slice_switch, FALSE);
	}
	if (in->dump_energy_slice_ypos>=in->ymeshpoints) in->dump_energy_slice_ypos=0;

	inp_search_int(sim,&inp,&(in->dump_energy_slice_zpos),"#dump_energy_slice_zpos");
	if (in->dump_energy_slice_xpos<0)
	{
			set_dump_status(sim,dump_energy_slice_switch, FALSE);
	}
	if (in->dump_energy_slice_zpos>=in->zmeshpoints) in->dump_energy_slice_zpos=0;


	set_dump_status(sim,dump_1d_slices, TRUE);

	inp_search_int(sim,&inp,&(in->dump_1d_slice_xpos),"#dump_1d_slice_xpos");
	if (in->dump_1d_slice_xpos<0) in->dump_1d_slice_xpos=-1;

	inp_search_int(sim,&inp,&(in->dump_1d_slice_zpos),"#dump_1d_slice_zpos");
	if (in->dump_1d_slice_zpos<0) in->dump_1d_slice_zpos=-1;


	dump=inp_search_english(sim,&inp,"#dump_verbose_electrical_solver_results");
	set_dump_status(sim,dump_1d_slices, dump);


	dump=inp_search_english(sim,&inp,"#dump_print_newtonerror");
	set_dump_status(sim,dump_print_newtonerror, dump);

	dump=inp_search_english(sim,&inp,"#dump_write_converge");
	set_dump_status(sim,dump_write_converge, dump);

	dump=inp_search_english(sim,&inp,"#dump_print_converge");
	set_dump_status(sim,dump_print_converge, dump);

	dump=inp_search_english(sim,&inp,"#dump_print_pos_error");
	set_dump_status(sim,dump_print_pos_error, dump);

	dump=inp_search_english(sim,&inp,"#dump_pl");
	set_dump_status(sim,dump_pl, dump);

	dump=inp_search_english(sim,&inp,"#dump_zip_files");
	set_dump_status(sim,dump_zip_files, dump);

	dump=inp_search_english(sim,&inp,"#dump_write_out_band_structure");
	set_dump_status(sim,dump_write_out_band_structure, dump);

	dump=inp_search_english(sim,&inp,"#dump_first_guess");
	set_dump_status(sim,dump_first_guess, dump);

	dump=inp_search_english(sim,&inp,"#dump_optical_probe");
	set_dump_status(sim,dump_optical_probe, dump);

	dump=inp_search_english(sim,&inp,"#dump_optical_probe_spectrum");
	set_dump_status(sim,dump_optical_probe_spectrum, dump);


	sim->log_level=inp_search_english(sim,&inp,"#dump_log_level");

	dump=inp_search_english(sim,&inp,"#dump_print_text");
	set_dump_status(sim,dump_print_text, dump);

	dump=inp_search_english(sim,&inp,"#dump_optics_summary");
	set_dump_status(sim,dump_optics_summary, dump);

	dump=inp_search_english(sim,&inp,"#dump_info_text");
	set_dump_status(sim,dump_info_text, dump);

	dump=inp_search_english(sim,&inp,"#dump_ray_trace_map");
	set_dump_status(sim,dump_ray_trace_map, dump);

	dump=inp_search_english(sim,&inp,"#dump_file_access_log");
	set_dump_status(sim,dump_file_access_log, dump);

	dump=inp_search_english(sim,&inp,"#dump_use_cache");
	set_dump_status(sim,dump_use_cache, dump);

	dump=inp_search_english(sim,&inp,"#dump_write_headers");
	set_dump_status(sim,dump_write_headers, dump);

	dump=inp_search_english(sim,&inp,"#dump_remove_dos_cache");
	set_dump_status(sim,dump_remove_dos_cache, dump);

	inp_free(sim,&inp);


}
