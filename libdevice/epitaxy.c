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

/** @file epitaxy.c
	@brief Load the epitaxy structure.
*/

#include <string.h>
#include "epitaxy.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include <cal_path.h>
#include <shape.h>
#include <contacts.h>

/**
 * @brief Load the epitaxy in from epitaxy.inp
 *
 */

void epitaxy_load_electrical_file(struct simulation *sim,struct epi_layer *layer)
{
	struct inp_file inp;
	char full_path[PATH_MAX];
	char file_name[100];

	strcpy(file_name,layer->electrical_file);
	strcat(file_name,".inp");

	join_path(2, full_path, get_input_path(sim), file_name);

	if ((inp_isfile(sim,full_path)==0)&&(strcmp(layer->electrical_file,"none")!=0))
	{

		inp_init(sim,&inp);

		inp_load(sim,&inp,full_path);
		inp_check(sim,&inp,1.0);

		inp_search_string(sim,&inp,layer->component,"#electrical_component");

		inp_search_gdouble(sim,&inp,&(layer->shunt),"#electrical_shunt");

		inp_search_gdouble(sim,&inp,&(layer->series),"#electrical_series");
		//printf("oh>>>>%d");
		//getchar();

		inp_search_gdouble(sim,&inp,&(layer->n_ideality),"#electrical_n");

		inp_search_gdouble(sim,&inp,&(layer->J0),"#electrical_J0");

		inp_free(sim,&inp);

	}else
	{
		layer->shunt=0.0;
		layer->series=0.0;
		layer->n_ideality=0.0;
		layer->J0=0.0;
	}
}

void epitaxy_load_pl_file(struct simulation *sim,char *pl_file, struct epi_layer *layer)
{
	struct inp_file inp;
	char full_path[PATH_MAX];
	char temp_path[PATH_MAX];
	char spectrum_file[100];
	char temp[100];
	int count=0;
	if (strcmp(pl_file,"none")!=0)
	{
		strcpy(temp_path,pl_file);

		strcat(temp_path,".inp");
		join_path(2, full_path, get_input_path(sim), temp_path);
		if (inp_isfile(sim,full_path)!=0)
		{
			ewe(sim,"pl file %s does not exist",full_path);
		}

		inp_init(sim,&inp);

		inp_load(sim,&inp,full_path);
		inp_check(sim,&inp,1.0);

		inp_search_gdouble(sim,&inp,&(layer->pl_fe_fh),"#pl_fe_fh");

		inp_search_gdouble(sim,&inp,&(layer->pl_fe_te),"#pl_fe_te");

		inp_search_gdouble(sim,&inp,&(layer->pl_te_fh),"#pl_te_fh");

		inp_search_gdouble(sim,&inp,&(layer->pl_th_fe),"#pl_th_fe");

		inp_search_gdouble(sim,&inp,&(layer->pl_fh_th),"#pl_fh_th");

		inp_search_string(sim,&inp,temp,"#pl_emission_enabled");

		layer->pl_enabled=english_to_bin(sim,temp);
		inp_search_string(sim,&inp,temp,"#pl_use_experimental_emission_spectra");
		layer->pl_use_experimental_emission_spectra=english_to_bin(sim,temp);

		inp_search_gdouble(sim,&inp,&(layer->pl_experimental_emission_efficiency),"#pl_experimental_emission_efficiency");


		inp_search_string(sim,&inp,spectrum_file,"#pl_input_spectrum");


		if (strcmp(spectrum_file,"none")==0)
		{

			strcpy(layer->pl_spectrum_file,"none");
			if (layer->pl_use_experimental_emission_spectra==TRUE)
			{
				ewe(sim,"You have not specified an emisison spectra");
			}
		}else
		{
			join_path(3, layer->pl_spectrum_file, sim->emission_path, spectrum_file,"spectra.inp");

			epitaxy_load_emission(sim,layer);

		}


		inp_free(sim,&inp);

	}else
	{
		layer->pl_fe_fh=0.0;
		layer->pl_fe_te=0.0;
		layer->pl_te_fh=0.0;
		layer->pl_th_fe=0.0;
		layer->pl_fh_th=0.0;
		layer->pl_enabled=FALSE;
		strcpy(layer->pl_spectrum_file,"none");
	}
	layer->photon_extract_eff=NULL;
	layer->avg_photon_extract_eff=0.0;

	strcpy(layer->pl_file,pl_file);

}

void epitaxy_load_dos_files(struct simulation *sim,struct epitaxy *in, char *dos_file,char *lumo_file,char *homo_file)
{
	char temp[20];
	char full_path[PATH_MAX];

	strcpy(temp,dos_file);
	strcat(temp,".inp");
	join_path(2, full_path, get_input_path(sim), temp);

	if (inp_isfile(sim,full_path)!=0)
	{
		ewe(sim,"dos file %s does not exist",full_path);
	}
	strcpy(in->layer[in->electrical_layers].dos_file,dos_file);


	if (strcmp(lumo_file,"none")!=0)
	{
		strcpy(temp,lumo_file);
		strcat(temp,".inp");
		join_path(2, full_path, get_input_path(sim), temp);
		if (inp_isfile(sim,full_path)!=0)
		{
			ewe(sim,"lumo file %s does not exist",full_path);
		}
	}
	strcpy(in->lumo_file[in->electrical_layers],lumo_file);

	if (strcmp(homo_file,"none")!=0)
	{
		strcpy(temp,homo_file);
		strcat(temp,".inp");
		join_path(2, full_path, get_input_path(sim), temp);
		if (inp_isfile(sim,full_path)!=0)
		{
			ewe(sim,"homo file %s does not exist",full_path);
		}
	}
	strcpy(in->homo_file[in->electrical_layers],homo_file);

	in->electrical_layers++;
}
void epitaxy_load(struct simulation *sim,struct epitaxy *in, char *file)
{
	int i;
	char temp[100];
	char dos_file[20];
	char electrical_file[20];
	char pl_file[20];
	char lumo_file[20];
	char homo_file[20];

	long double y_pos=0.0;

	struct inp_file inp;
	in->electrical_layers=0;
	in->device_start=-1.0;
	in->device_stop=-1.0;
	char *test;
	inp_init(sim,&inp);
	if (inp_load(sim, &inp , file)!=0)
	{
		ewe(sim,"Epitaxy: I can't find file %s\n",file);
	}

	inp_check(sim,&inp,1.41);
	inp_reset_read(sim,&inp);
	test=inp_get_string(sim,&inp);

	sscanf(inp_get_string(sim,&inp),"%d",&(in->layers));

	if (in->layers>20)
	{
		ewe(sim,"Too many material layers\n");
	}

	if (in->layers<1)
	{
		ewe(sim,"No material layers\n");
	}

	for (i=0;i<in->layers;i++)
	{
		in->layer[i].layer_number=i;

		//name
		test=inp_get_string(sim,&inp);
		strcpy(in->layer[i].name,inp_get_string(sim,&inp));

		//type	(not yet used)
		test=inp_get_string(sim,&inp);
		strcpy(temp,inp_get_string(sim,&inp));
		in->layer[i].layer_type=english_to_bin(sim,temp);

		//width
		inp_get_string(sim,&inp);
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->layer[i].width));
		in->layer[i].width=fabs(in->layer[i].width);

		inp_get_string(sim,&inp);	//token
		strcpy(in->mat_file[i],inp_get_string(sim,&inp));
		assert_platform_path(in->mat_file[i]);

		inp_get_string(sim,&inp);	//token
		strcpy(dos_file,inp_get_string(sim,&inp));
		//strcpy(in->layer[i].layer_type,dos_file);


		inp_get_string(sim,&inp);	//token
		strcpy(pl_file,inp_get_string(sim,&inp));

		epitaxy_load_pl_file(sim,pl_file,&(in->layer[i]));

		inp_get_string(sim,&inp);	//token
		strcpy(in->shape_file[i],inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//token
		strcpy(lumo_file,inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//token
		strcpy(homo_file,inp_get_string(sim,&inp));

		//electrical file
		test=inp_get_string(sim,&inp);
		strcpy(in->layer[i].electrical_file,inp_get_string(sim,&inp));

		epitaxy_load_electrical_file(sim,&(in->layer[i]));

		char temp[20];
		char full_path[PATH_MAX];
		if (strcmp_begin(dos_file,"dos")==0)
		{
			if (in->device_start==-1.0)
			{
				in->device_start=y_pos;
			}

			in->layer[i].electrical_layer=TRUE;
			epitaxy_load_dos_files(sim,in,dos_file,lumo_file,homo_file);
			sprintf(electrical_file,"electrical%s.inp",(dos_file+3));
		}else
		{
			in->layer[i].electrical_layer=FALSE;
		}

		in->layer[i].y_start=y_pos;
		y_pos+=in->layer[i].width;
		in->layer[i].y_stop=y_pos;

	}

	in->device_stop=epitaxy_get_device_stop(in);


	char * ver = inp_get_string(sim,&inp);
	if (strcmp(ver,"#ver")!=0)
	{
			ewe(sim,"No #ver tag found in file\n");
	}

	inp_free(sim,&inp);

	epitaxy_load_materials(sim,in);


	epitaxy_shapes_load(sim,in);
}

/**
 * @brief Get the height of the layers which are electrically active.
 *
 */
gdouble epitaxy_get_electrical_length(struct epitaxy *in)
{
int i=0;
gdouble tot=0.0;

for (i=0;i<in->layers;i++)
{
	if (in->layer[i].electrical_layer==TRUE)
	{
		tot+=in->layer[i].width;
	}
}
//if (tot>300e-9)
//{
//	ewe(sim,"Can't simulate structures bigger than 300 nm\n");
//}
return tot;
}

/**
 * @brief Get the height of the layers which are optically active. i.e. the full height of the device.
 *
 */
gdouble epitaxy_get_optical_length(struct epitaxy *in)
{
int i=0;
gdouble tot=0.0;

for (i=0;i<in->layers;i++)
{
	tot+=in->layer[i].width;
}

return tot;
}

/**
 * @brief Return the index of a layer in the device.
 *
 */
int epitaxy_get_layer(struct epitaxy *in,gdouble pos)
{
int i=0;
long double layer_end=0.0;//-in->device_start;
for (i=0;i<in->layers;i++)
{
	layer_end+=in->layer[i].width;

	if (pos<=layer_end)
	{
		return i;
	}

}

return -1;
}


/**
 * @brief Return the electrical layer number of a layer from a y cordinate.
 *
 */
int epitaxy_get_electrical_material_layer(struct epitaxy *in,gdouble pos)
{
int i=0;
gdouble layer_end=0.0;
int electrical_layer=0;
for (i=0;i<in->layers;i++)
{
	if (in->layer[i].electrical_layer==TRUE)
	{
		layer_end+=in->layer[i].width;

		if (pos<layer_end)
		{
			return electrical_layer;
		}
		electrical_layer++;
	}

}

return -1;
}

/**
 * @brief Return the layer number of a layer from a y cordinate.
 *
 */
int epitaxy_get_epitaxy_layer_using_electrical_pos(struct epitaxy *in,gdouble pos)
{
int i=0;
gdouble layer_end=0.0;
int electrical_layer=0;
for (i=0;i<in->layers;i++)
{
	if (in->layer[i].electrical_layer==TRUE)
	{
		layer_end+=in->layer[i].width;

		if (pos<layer_end)
		{
			return i;
		}
		electrical_layer++;
	}

}

return -1;
}

/**
 * @brief Return the position of the start of the first electrical layer in the device.
 *
 */
gdouble epitaxy_get_device_start(struct epitaxy *in)
{
int i=0;
gdouble pos=0.0;
for (i=0;i<in->layers;i++)
{

	if (in->layer[i].electrical_layer==TRUE)
	{
		return pos;
	}
	pos+=in->layer[i].width;

}

return -1;
}

/**
 * @brief Return the position of the end of the electrical layers in the device.
 *
 */
long double epitaxy_get_device_stop(struct epitaxy *in)
{
int i=0;
gdouble pos=0.0;
int found=FALSE;
for (i=0;i<in->layers;i++)
{

	if (in->layer[i].electrical_layer==TRUE)
	{
		found=TRUE;
	}

	if ((in->layer[i].electrical_layer==FALSE)&&(found==TRUE))
	{
		return pos;
	}

	pos+=in->layer[i].width;
}

if (found==TRUE)
{
	return pos;
}


return -1;
}

/**
 * @brief Return the index of the layer where the device starts.
 *
 */
gdouble epitaxy_get_device_start_i(struct epitaxy *in)
{
int i=0;
for (i=0;i<in->layers;i++)
{

	if (in->layer[i].electrical_layer==TRUE)
	{
		return i;
	}

}

return -1;
}

/**
 * @brief Free the epitaxy structure.
 *
 */
void epitaxy_free(struct simulation *sim,struct epitaxy *in)
{
	epitaxy_free_materials(in);
	shape_free_materials(sim,in);
}

void epitaxy_mask(struct simulation *sim,struct device *dev)
{
int x=0;
int y=0;
int z=0;
int c=0;
int l=0;
struct dimensions *dim=&(dev->ns.dim);
struct epitaxy *epi=&(dev->my_epitaxy);

if ((dim->xlen==1)&&(dim->zlen==1))
{
	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				dev->mask[z][x][y]=TRUE;
			}
		}
	}
	return;
}

	for (z=0;z<dim->zlen;z++)
	{
		for (x=0;x<dim->xlen;x++)
		{
			for (y=0;y<dim->ylen;y++)
			{
				dev->mask[z][x][y]=FALSE;

				l=dev->imat_epitaxy[z][x][y];
				if (l==0)	//first layer
				{
					if (dev->my_epitaxy.layer[l].layer_type==LAYER_CONTACT)
					{
						for (c=0;c<dev->ncontacts;c++)
						{
							if (dev->contacts[c].position==TOP)
							{
								if (contact_within_zx(&(dev->contacts[c]), dim->zmesh[z], dim->xmesh[x])==0)
								{
									dev->mask[z][x][y]=TRUE;
								}
							}
						}
					}
				}else
				{
					dev->mask[z][x][y]=TRUE;
				}
			}
		}
	}

return;
}

