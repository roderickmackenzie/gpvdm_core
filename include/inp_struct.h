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

/** @file inp_struct.h
@brief structure to hold input files.
*/


#ifndef inpstruct_h
#define inpstruct_h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "advmath.h"

struct inp_file
{
char *data;
long fsize;
char full_name[1000];
int pos;
int edited;
};

struct inp_list
{
	char **names;
	int len;
};
#endif
