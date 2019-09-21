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

/** @file basic_math.c
	@brief Math functions
*/

#include <math.h>


long double log_delta(long double a_in,long double b_in)
{
	long double a=fabsl(a_in);
	long double b=fabsl(b_in);

	long double ratio=0.0;

	if (a<b)
	{
		b=a;
		a=fabsl(b_in);
	}


	if (b==a)
	{
		return 0.0;
	}

	if ((b==0)&&(a==0))
	{
		return 0;
	}

	if (b==0)
	{
		return fabsl(log10l(a));
	}
	//printf("one %Le %Le %Le\n",a,b,fabsl(log10l(a/b)));
	ratio=fabsl(log10l(a/b));
	
	return ratio;
}

