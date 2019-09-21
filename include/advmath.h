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

/** @file advmath.h
	@brief Math functions which have been moved away from the main code for compiler optimization.
*/
#ifndef advmath_h
#define advmath_h
#include <math.h>

#define gdouble long double
#define gpow powl
#define gcabs cabsl
#define gcreal creall
#define gcimag cimagl
#define gfabs fabsl
#define gexp expl
#define gsqrt sqrtl

long double B(long double x);
long double dB(long double x);
long double log_delta(long double a_in,long double b_in);
#endif
