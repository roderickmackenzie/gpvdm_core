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

/** @file circuit.h
@brief Header files for nodal analysis
*/

#ifndef circuit_struct_h
#define circuit_struct_h
#include <sim_struct.h>

#define CIR_KNOWN 0
#define CIR_UNKNOWN 1
#define CIR_CHANGE_X 2

struct circuit_config_line
{
	char component[100];
	int x0;
	int y0;
	int z0;
	int x1;
	int y1;
	int z1;
	double R;
	double C;
	double L;
	double nid;
	double J0;
};

struct circuit_link
{
	int start;
	int stop;
	char type;
	double R;
	double L;
	double C;
	double n0;
	double J0;
	double Jsc;
	double j;

	//only used for diodes and to figure out where they are.
	int layer;
	int x;
	int z;
	int id;

};

struct circuit_node
{
	double V;
	double V_last;
	int type;
	int matrix_pos;
	int z;
	int x;
	int y;

	double z_pos;
	double x_pos;
	double y_pos;

	int links[10];
	int nlinks;
};

struct circuit
{
	int nodes_len;
	int links_len;
	int nodes_max;
	int links_max;
	int unknowns;
	struct matrix mx;
	struct circuit_node *nodes;
	struct circuit_link *links;
	int config_nlines;
	struct circuit_config_line * config_lines;
	int abstract_circuit;		//The circuit does not follow a device structure
	int enabled;
};

#endif
