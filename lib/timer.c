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

/** @file timer.c
	@brief Timer.
*/


#include <sys/time.h>
#include <stdio.h>

static double start_times[10];

void timer_init(int timer)
{
struct timeval result;
gettimeofday (&result, NULL);
start_times[timer]=result.tv_sec + result.tv_usec/1000000.0;
}

double timer_get_time(int timer)
{
struct timeval result;
gettimeofday (&result, NULL);
double cur_time=result.tv_sec + result.tv_usec/1000000.0;
return cur_time-start_times[timer];
}

