//
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/** @file memory_basic.c
@brief memory functions for 3D arrays
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lang.h>
#include "sim.h"
#include "dump.h"
#include "mesh.h"
#include <math.h>
#include "log.h"
#include <solver_interface.h>
#include "memory.h"


void malloc_light_zxyl_long_double_complex(struct dim_light *dim, long double complex * (****var))
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;

	*var = (long double complex ****) malloc(dim->zlen * sizeof(long double complex ***));

	for (z = 0; z < dim->zlen; z++)
	{
		(*var)[z] = (long double complex ***) malloc(dim->xlen * sizeof(long double complex**));
		for (x = 0; x < dim->xlen; x++)
		{
			(*var)[z][x] = (long double complex **) malloc(dim->ylen * sizeof(long double complex*));
			for (y = 0; y < dim->ylen; y++)
			{
				(*var)[z][x][y] = (long double complex *) malloc(dim->llen * sizeof(long double complex));
				memset((*var)[z][x][y], 0, dim->llen * sizeof(long double complex));

			}
		}
	}

}



void free_light_zxyl_long_double_complex(struct dim_light *dim, long double complex * (****in_var))
{
	int x=0;
	int y=0;
	int z=0;
	int l=0;

	long double complex ****var=*in_var;
	if (var==NULL)
	{
		return;
	}

	for (z = 0; z < dim->zlen; z++)
	{

		for (x = 0; x < dim->xlen; x++)
		{
			for (y = 0; y < dim->ylen; y++)
			{
				free(var[z][x][y]);
			}
			free(var[z][x]);
		}
		free(var[z]);
	}

	free(var);

	*in_var=NULL;

}

