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

#include <stdio.h>
#include <ray.h>
#include <gpvdm_const.h>
#include <math.h>
#include <stdlib.h>
#include <cal_path.h>
#include <log.h>

/** @file wavelength_to_rgb.c
	@brief Turn a wavelegngth to an RGB color.
*/

void wavelength_to_rgb(int *r,int *g,int *b,double wavelength)
{
	double gamma=0.80;
	double intensity_max=1.0;
	double factor=0.0;

	double red=0.0;
	double green=0.0;
	double blue=0.0;

	if ((wavelength >= 380e-9) && (wavelength<440e-9))
	{
		red = -(wavelength - 440e-9) / (440e-9 - 380e-9);
		green = 0.0;
		blue = 1.0;
	}else
	if ((wavelength >= 440e-9) &&  (wavelength<490e-9))
	{
		red = 0.0;
		green = (wavelength - 440e-9) / (490e-9 - 440e-9);
		blue = 1.0;
	}else
	if ((wavelength >= 490e-9) &&  (wavelength<510e-9))
	{
		red = 0.0;
		green = 1.0;
		blue = -(wavelength - 510e-9) / (510e-9 - 490e-9);
	}else
	if ((wavelength >= 510e-9) &&  (wavelength<580e-9))
	{
		red = (wavelength - 510e-9) / (580e-9 - 510e-9);
		green = 1.0;
		blue = 0.0;
	}else
	if ((wavelength >= 580e-9) &&  (wavelength<645e-9))
	{
		red = 1.0;
		green = -(wavelength - 645e-9) / (645e-9 - 580e-9);
		blue = 0.0;
	}else
	if ((wavelength >= 645e-9) &&  (wavelength<781e-9))
	{
		red = 1.0;
		green = 0.0;
		blue = 0.0;
	}else
	{
		red = 0.0;
		green = 0.0;
		blue = 0.0;
	}

	if (wavelength<420e-9)
	{
		factor = 0.3 + 0.7*(wavelength - 380e-9) / (420e-9 - 380e-9);
	}else
	if ((wavelength >= 420e-9) && (wavelength<701e-9))
	{
		factor = 1.0;
	}else
	if ((wavelength >= 701e-9) &&  (wavelength<781e-9))
	{
		factor = 0.3 + 0.7*(780e-9 - wavelength) / (780e-9 - 700e-9);
	}else
	{
		factor = 0.0;
	}

	if (red != 0.0)
	{
		red = (intensity_max * pow(red * factor, gamma));
	}

	if (green != 0)
	{
		green = (intensity_max * pow(green * factor, gamma));
	}

	if (blue != 0)
	{
		blue = (intensity_max * pow(blue * factor, gamma));
	}


	if ((red==0.0)&&(green==0.0)&&(blue==0.0))
	{
		red=0.392634;
		green=0.000000;
		blue=0.397315;
	}

	*r=255.0*red;
	*g=255.0*green;
	*b=255.0*blue;
}
