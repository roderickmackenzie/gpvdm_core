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

/** @file checksum.h
	@brief MD5 type checksums.
*/


#ifndef md5_h
#define md5_h
#include <stdint.h>
#include <sim_struct.h>

struct md5
{
	uint32_t a0;
	uint32_t b0;
	uint32_t c0;
	uint32_t d0;
};

uint32_t leftrotate (uint32_t x, uint32_t c);
void md5_init(struct md5* in);
void md5_update(struct md5* in,char *data,int len);
void md5_to_str(char *out,struct md5 *in);

#endif
