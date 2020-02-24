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

/** @file contacts.c
@brief backend to handle complex contacts
*/

#include <string.h>
#include "epitaxy.h"
#include "inp.h"
#include "util.h"
#include "gpvdm_const.h"
#include <cal_path.h>
#include "contacts.h"
#include <log.h>
#include <dump.h>
#include <shape.h>
#include <shape_struct.h>

int contacts_get_lcharge_type(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].position==TOP)
		{
			return in->contacts[i].charge_type;
		}
	}

ewe(sim,"No top contact found\n");
}

int contacts_get_rcharge_type(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].position==BOTTOM)
		{
			return in->contacts[i].charge_type;
		}
	}

ewe(sim,"No bottom contact found\n");
}

long double contacts_get_lcharge(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].position==TOP)
		{
			return in->contacts[i].np;
		}
	}

ewe(sim,"No top contact found\n");
}

long double contacts_get_rcharge(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].position==BOTTOM)
		{
			return in->contacts[i].np;
		}
	}
ewe(sim,"No bottom contact found\n");
}

void contacts_time_step(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		in->contacts[i].voltage_last=in->contacts[i].voltage;
	}
}

int contacts_itterate_to_desired_voltage(struct simulation *sim,struct device *in,char *text)
{
int i;
char temp[400];
static long double dV=0.1;
int up=TRUE;
int changed=FALSE;
struct newton_state *ns=&(in->ns);

	strcpy(text,"Ramping voltage: ");

	for (i=0;i<in->ncontacts;i++)
	{

		up=TRUE;
		if (in->contacts[i].voltage_want!=in->contacts[i].voltage)
		{
			changed=TRUE;
			if ((in->contacts[i].voltage_want-in->contacts[i].voltage)<0.0)
			{
				up=FALSE;
			}

			if (ns->last_ittr<16)
			{
				dV=dV*(1.0+0.001);
			}

			if (ns->last_ittr>18)
			{
				dV*=(1.0-0.001);
			}

			if (dV<0.1)
			{
				dV=0.1;
			}

			/*if (in->contacts[i].voltage<5.0)
			{
				dV=0.1;
			}else
			{
				dV=0.1+fabs(in->contacts[i].voltage)*0.005;
			}*/

			if (up==TRUE)
			{
				if ((in->contacts[i].voltage+dV)>=in->contacts[i].voltage_want)
				{
					in->contacts[i].voltage=in->contacts[i].voltage_want;
				}else
				{
					in->contacts[i].voltage+=dV;
				}

			}else
			{
				if ((in->contacts[i].voltage-dV)<=in->contacts[i].voltage_want)
				{
					in->contacts[i].voltage=in->contacts[i].voltage_want;
				}else
				{
					in->contacts[i].voltage-=dV;
				}
			}

			printf_log(sim,"Ramping: %s %.2Lf %.2Lf dV=%Lf ittr=%d\n",in->contacts[i].name,in->contacts[i].voltage,in->contacts[i].voltage_want,dV,ns->last_ittr);

			sprintf(temp,"%s %.2Lf V/%.2Lf V",in->contacts[i].name,in->contacts[i].voltage,in->contacts[i].voltage_want);
			strcat(text,temp);
		}


	}

	if (changed==TRUE)
	{
		contacts_update(sim,in);
	}

	return changed;
}

void contacts_load(struct simulation *sim,struct device *in)
{
	int i;
	struct inp_file inp;
	in->ncontacts=0;

	inp_init(sim,&inp);
	if (inp_load(sim, &inp , "contacts.inp")!=0)
	{
		ewe(sim,"Can't open the file contacts\n");
	}

	inp_check(sim,&inp,1.3);
	inp_reset_read(sim,&inp);
	inp_get_string(sim,&inp);
	sscanf(inp_get_string(sim,&inp),"%d",&(in->ncontacts));

	if (in->ncontacts>10)
	{
		ewe(sim,"Too many contacts\n");
	}

	int active=FALSE;
	long double ingress=0.0;

	struct shape *s;

	for (i=0;i<in->ncontacts;i++)
	{

		inp_get_string(sim,&inp);	//name
		strcpy(in->contacts[i].name,inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//position
		in->contacts[i].position=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//active contact
		in->contacts[i].active=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//voltage
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].voltage_want));
		in->contacts[i].voltage=0.0;
		in->contacts[i].voltage_last=in->contacts[i].voltage;
		//printf("%Le\n",in->contacts[i].voltage_want);
		//getchar();

		inp_get_string(sim,&inp);	//np
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].np));
		in->contacts[i].np=fabs(in->contacts[i].np);

		inp_get_string(sim,&inp);	//charge_type
		in->contacts[i].charge_type=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//shape_file_name
		strcpy(in->contacts[i].shape_file_name,inp_get_string(sim,&inp));

		s=shape_load_file(sim,&(in->my_epitaxy),&(in->contacts[i].shape),in->contacts[i].shape_file_name);
		s->nx=1;
		s->nz=1;
		s->dz=in->zlen;
		strcpy(s->name,in->contacts[i].name);

		inp_get_string(sim,&inp);	//resistance sq
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].Rcontact));
		in->contacts[i].Rcontact=fabs(in->contacts[i].Rcontact);

		inp_get_string(sim,&inp);	//ingress
		sscanf(inp_get_string(sim,&inp),"%Le",&(ingress));
		ingress=fabs(ingress);

		inp_get_string(sim,&inp);	//contact type
		in->contacts[i].type=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//ve0
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].ve0));
		in->contacts[i].ve0=fabs(in->contacts[i].ve0);

		inp_get_string(sim,&inp);	//ve0
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].vh0));
		in->contacts[i].vh0=fabs(in->contacts[i].vh0);

		if (in->contacts[i].position==LEFT)
		{
			s->dx=ingress;
		}
	}

	char * ver = inp_get_string(sim,&inp);
	if (strcmp(ver,"#ver")!=0)
	{
			ewe(sim,"No #ver tag found in file\n");
	}

	inp_free(sim,&inp);

	contacts_update(sim,in);
	contact_set_flip_current(sim,in);

	contacts_cal_area(sim,in);

}

void contacts_force_to_zero(struct simulation *sim,struct device *in)
{
int y;
int x;
int z;
struct dimensions *dim=&in->ns.dim;

for (z=0;z<dim->zlen;z++)
{
	for (x=0;x<dim->xlen;x++)
	{
			in->Vapplied_y0[z][x]=0.0;
			in->Vapplied_y1[z][x]=0.0;
	}

	for (y=0;y<dim->ylen;y++)
	{
			in->Vapplied_x0[z][y]=0.0;
			in->Vapplied_x1[z][y]=0.0;
	}

}

}


void contacts_dump(struct simulation *sim,struct device *in)
{
	int i;
	int y=0;
	int x=0;
	int z=0;

	struct dimensions *dim=&in->ns.dim;

	if (get_dump_status(sim,dump_print_text)==TRUE)
	{
		for (i=0;i<in->ncontacts;i++)
		{
			printf_log(sim,"%-10s\tV=%Le\tA=%Le\n",in->contacts[i].name,in->contacts[i].voltage,in->contacts[i].area);
		}

	}



	printf("top-btm\n");
	for (z = 0; z < dim->zlen; z++)
	{
		for (x = 0; x < dim->xlen; x++)
		{
			printf("%d %d\n",in->n_contact_y0[z][x],in->n_contact_y1[z][x]);

			//printf("%d %.2Le %.2Lf %.2Lf (%.2Lf %.2Lf %Le %Le %Le %Le)\n",x,dim->xmesh[x],in->Vapplied_y0[z][x],in->Vapplied_y1[z][x],in->V_y0[z][x],in->V_y1[z][x],in->electrons_y0[z][x],in->holes_y0[z][x],in->electrons_y1[z][x],in->holes_y1[z][x]);
		}
	}

	printf("left-right\n");
	for (z = 0; z < dim->zlen; z++)
	{
		for (y = 0; y < dim->ylen; y++)
		{
			printf("%d %d\n",in->n_contact_x0[z][y],in->n_contact_x1[z][y]);
			//printf("%d %.2Le %.2Lf %.2Lf (%.2Lf %.2Lf %Le %Le %Le %Le)\n",y,dim->ymesh[y],in->Vapplied_x0[z][y],in->Vapplied_x1[z][y],in->V_x0[z][y],in->V_x1[z][y],in->electrons_x0[z][y],in->holes_x0[z][y],in->electrons_x1[z][y],in->holes_x1[z][y]);
		}
	}
}

int contact_within_zx(struct contact *c, long double z,long double x)
{
	if (x>=c->shape.x0)
	{
		if (x<(c->shape.x0+c->shape.dx))
		{
			return 0;
		}
	}

	return -1;
}

int contact_within_zy(struct contact *c, long double z, long double y)
{

	if (y>=c->shape.y0)
	{
		if (y<(c->shape.y0+c->shape.dy))
		{
			return 0;
		}
	}

	return -1;
}

void contacts_update(struct simulation *sim,struct device *in)
{
int i;
int y;
int x;
int z;
int n;
int found=FALSE;

gdouble value=0.0;
struct newton_state *ns=&in->ns;
struct dimensions *dim=&in->ns.dim;
struct contact *c;

if (in->Vapplied_y0==NULL) return;

if (dim->xlen==1)
{
	for (z=0;z<dim->zlen;z++)
	{

		in->passivate_y0[z][0]=FALSE;
		in->passivate_y1[z][0]=FALSE;

		for (i=0;i<in->ncontacts;i++)
		{
			c=&(in->contacts[i]);
			if ((c->position==TOP)&&(c->active==TRUE))
			{
				in->Vapplied_y0[z][0]=c->voltage;
				in->n_contact_y0[z][0]=i;
			}else
			{
				in->Vapplied_y1[z][0]=c->voltage;
				in->n_contact_y1[z][0]=i;
			}
		}
	}

	//contacts_dump(sim,in);
	//getchar();
	return;
}

//Reset contacts
for (z=0;z<dim->zlen;z++)
{
	for (x=0;x<dim->xlen;x++)
	{
		in->Vapplied_y0[z][x]=0.0;
		in->Vapplied_y1[z][x]=0.0;

		in->n_contact_y0[z][x]=-1;
		in->n_contact_y1[z][x]=-1;

		in->passivate_y0[z][x]=TRUE;
		in->passivate_y1[z][x]=TRUE;
	}

	for (y=0;y<dim->ylen;y++)
	{
		in->Vapplied_x0[z][y]=0.0;
		in->Vapplied_x1[z][y]=0.0;

		in->n_contact_x0[z][y]=-1;
		in->n_contact_x1[z][y]=-1;

		in->passivate_x0[z][y]=TRUE;
		in->passivate_x1[z][y]=TRUE;
	}

}

	for (i=0;i<in->ncontacts;i++)
	{
		c=&(in->contacts[i]);
		if (c->position==TOP)
		{
			for (z=0;z<dim->zlen;z++)
			{
				for (x=0;x<dim->xlen;x++)
				{

					if (contact_within_zx(c, dim->zmesh[z],dim->xmesh[x])==0)
					{

						if (in->n_contact_y0[z][x]!=-1)
						{
							ewe(sim,"You have overlapping contacts\n");
						}

						in->Vapplied_y0[z][x]=c->voltage;
						in->n_contact_y0[z][x]=i;
						in->passivate_y0[z][x]=FALSE;
					}
				}
			}
		}else
		if (c->position==BOTTOM)
		{
			//printf("btm %d\n",in->contacts[i].position);
			for (z=0;z<dim->zlen;z++)
			{
				for (x=0;x<dim->xlen;x++)
				{

					if (contact_within_zx(c, dim->zmesh[z],dim->xmesh[x])==0)
					{

						if (in->n_contact_y1[z][x]!=-1)
						{
							ewe(sim,"You have overlapping contacts\n");
						}

						in->Vapplied_y1[z][x]=c->voltage;
						in->n_contact_y1[z][x]=i;
						in->passivate_y1[z][x]=FALSE;
					}
				}
			}
		}else
		if (c->position==LEFT)
		{
			for (z=0;z<dim->zlen;z++)
			{
				for (y=0;y<dim->ylen;y++)
				{


					if (contact_within_zy(c, dim->zmesh[z], dim->ymesh[y]+in->my_epitaxy.device_start)==0)
					{

						if (in->n_contact_x0[z][y]!=-1)
						{
							ewe(sim,"You have overlapping contacts\n");
						}

						in->Vapplied_x0[z][y]=c->voltage;
						in->n_contact_x0[z][y]=i;
						in->passivate_x0[z][y]=FALSE;
					}
				}
			}
		}else
		if (c->position==RIGHT)
		{
			for (z=0;z<dim->zlen;z++)
			{
				for (y=0;y<dim->ylen;y++)
				{

					if (contact_within_zy(c, dim->zmesh[z], dim->ymesh[y]+in->my_epitaxy.device_start)==0)
					{

						if (in->n_contact_x1[z][y]!=-1)
						{
							ewe(sim,"You have overlapping contacts\n");
						}

						in->Vapplied_x1[z][y]=c->voltage;
						in->n_contact_x1[z][y]=i;
						in->passivate_x1[z][y]=FALSE;
					}
				}
			}
		}
	}

//contacts_dump(sim,in);
//getchar();
}

gdouble contact_get_voltage_last(struct simulation *sim,struct device *in,int contact)
{
	return in->contacts[contact].voltage_last;
}

gdouble contact_get_voltage(struct simulation *sim,struct device *in,int contact)
{
	return in->contacts[contact].voltage;
}

void contact_set_voltage(struct simulation *sim,struct device *in,int contact,gdouble voltage)
{
	in->contacts[contact].voltage=voltage;
	contacts_update(sim,in);
}

void contact_set_wanted_active_contact_voltage(struct simulation *sim,struct device *in,gdouble voltage)
{
	int i=0;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			in->contacts[i].voltage_want=voltage;
		}
	}

	contacts_update(sim,in);

}

void contact_set_active_contact_voltage(struct simulation *sim,struct device *in,gdouble voltage)
{
	int i=0;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			in->contacts[i].voltage=voltage;
		}
	}

	contacts_update(sim,in);
	//contacts_dump(sim,in);
}

void contact_set_flip_current(struct simulation *sim,struct device *in)
{
	int i=0;
	in->flip_current=1.0;
	if (in->ncontacts==2)
	{
		for (i=0;i<in->ncontacts;i++)
		{
			//printf("%d %d %d\n",in->contacts[i].active,in->contacts[i].position,BOTTOM);
			//getchar();

			if ((in->contacts[i].active==TRUE)&&(in->contacts[i].position==BOTTOM))
			{
				in->flip_current=-1.0;
			}
		}
	}
}

gdouble contact_get_active_contact_voltage_last(struct simulation *sim,struct device *in)
{
	int i=0;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			return in->contacts[i].voltage_last;
		}
	}
}

int contact_get_active_contact_index(struct simulation *sim,struct device *in)
{
	int i=0;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			return i;
		}
	}


}


long double contact_get_active_contact_voltage(struct simulation *sim,struct device *in)
{
	int i=0;

	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			return in->contacts[i].voltage;
		}
	}


}

void contact_set_all_voltages(struct simulation *sim,struct device *in,gdouble voltage)
{
int i;
	for (i=0;i<in->ncontacts;i++)
	{
		in->contacts[i].voltage=voltage;
	}

	contacts_update(sim,in);
}

long double contacts_get_Jleft(struct device *in)
{
int i;
int x;
int z;

long double tot=0.0;
long double count=0.0;
struct dimensions *dim=&in->ns.dim;

for (x=0;x<dim->xlen;x++)
{
		for (z=0;z<dim->zlen;z++)
		{
			if (in->n_contact_y0[z][x]>=0)
			{
				tot+=in->Jp_y0[z][x]+in->Jn_y0[z][x];
				count=count+1.0;						//this will need updating for meshes which change
			}
		}
}

tot=tot/count;
tot*=in->flip_current;

return tot;
}

long double contacts_get_Jright(struct device *in)
{
int i;
int x;
int z;

long double tot=0.0;
long double count=0.0;
struct dimensions *dim=&in->ns.dim;

for (x=0;x<dim->xlen;x++)
{
		for (z=0;z<dim->zlen;z++)
		{
			if (in->n_contact_y1[z][x]>=0)
			{
				tot+=in->Jp_y1[z][x]+in->Jn_y1[z][x];
				count=count+1.0;
			}
		}
}

tot=tot/count;

tot*=in->flip_current;

return tot;
}

int contacts_get_active_contact_left_right(struct device *in)
{
int i;
	for (i=0;i<in->ncontacts;i++)
	{
		if (in->contacts[i].active==TRUE)
		{
			if (in->contacts[i].position==TOP)
			{
				return LEFT;
			}else
			{
				return RIGHT;
			}
		}

	}

return -1;
}

void contacts_cal_area(struct simulation *sim,struct device *in)
{
int i;
int x;
int z;

if (in->n_contact_y1==NULL) return;

struct dimensions *dim=&in->ns.dim;

for (i=0;i<in->ncontacts;i++)
{
	in->contacts[i].area=0.0;
}

for (x=0;x<dim->xlen;x++)
{
		for (z=0;z<dim->zlen;z++)
		{
			i=in->n_contact_y1[z][x];
			if (i!=-1)
			{
				in->contacts[i].area+=dim->dx[x]*dim->dz[z];
			}

			i=in->n_contact_y0[z][x];
			if (i!=-1)
			{
				in->contacts[i].area+=dim->dx[x]*dim->dz[z];
			}

		}
}

}

void contacts_detailed_dump(struct device *in)
{
int i;
int x;
int z;
struct dimensions *dim=&in->ns.dim;

for (x=0;x<dim->xlen;x++)
{
		for (z=0;z<dim->zlen;z++)
		{
			printf("%d %Le %Le\n",in->n_contact_y0[z][x],in->Jp_y0[z][x],in->Jn_y0[z][x]);
		}
}

}

//Average the current over both contacts
long double contacts_get_J(struct device *in, int n)
{
int i;
int x;
int z;

long double tot=0.0;
long double count=0.0;
struct dimensions *dim=&in->ns.dim;


for (x=0;x<dim->xlen;x++)
{
		for (z=0;z<dim->zlen;z++)
		{
			for (i=0;i<in->ncontacts;i++)
			{
				if (in->n_contact_y1[z][x]==n)
				{
					tot+=in->Jp_y1[z][x]+in->Jn_y1[z][x];
					count=count+1.0;						//this will need updating for meshes which change
				}

				if (in->n_contact_y0[z][x]==n)
				{
					tot+=in->Jp_y0[z][x]+in->Jn_y0[z][x];
					count=count+1.0;						//this will need updating for meshes which change
				}
			}
		}
}

tot=tot/count;

tot*=in->flip_current;

return tot;
}

void contacts_passivate(struct simulation *sim,struct device *in)
{
int i;
int x;
int y;
int z;
return;
struct newton_state *ns=&in->ns;
struct dimensions *dim=&in->ns.dim;

//passivate under each contact
for (x=0;x<dim->xlen;x++)
{
	for (y=0;y<dim->ylen;y++)
	{
		for (z=0;z<dim->zlen;z++)
		{

			for (i=0;i<in->ncontacts;i++)
			{
				if (in->contacts[i].position==TOP)
				{
					//if ((in->ylen-dim->ymesh[y]<=in->contacts[i].depth)&&(dim->xmesh[x]>in->contacts[i].start)&&(dim->xmesh[x]<in->contacts[i].start+in->contacts[i].width))
					{
						in->mun[z][x][y]=1e-15;
						in->mup[z][x][y]=1e-15;
					}
				}

				if (in->contacts[i].position==BOTTOM)
				{
					//if ((dim->ymesh[y]<=in->contacts[i].depth)&&(dim->xmesh[x]>in->contacts[i].start)&&(dim->xmesh[x]<in->contacts[i].start+in->contacts[i].width))
					{
						in->mun[z][x][y]=1e-15;
						in->mup[z][x][y]=1e-15;
					}
				}
			}
		}
	}
}

for (x=0;x<dim->xlen;x++)
{
	for (z=0;z<dim->zlen;z++)
	{
		i=in->n_contact_y1[z][x];
		if (i==-1)
		{
			in->mun[z][x][dim->ylen-1]=1e-15;
			in->mup[z][x][dim->ylen-1]=1e-15;
		}

		i=in->n_contact_y0[z][x];
		if (i==-1)
		{
			in->mun[z][x][0]=1e-15;
			in->mup[z][x][0]=1e-15;
		}
	}
}

}

void contacts_free(struct simulation *sim,struct device *in)
{
	int i;

	for (i=0;i<in->ncontacts;i++)
	{
		shape_free(sim,&(in->contacts[i].shape));
	}
}

void contacts_ingress(struct simulation *sim,struct device *in)
{
	int x=0;
	int y=0;
	int z=0;
	int c=0;
	struct newton_state *ns=&(in->ns);
	struct dimensions *dim=&in->ns.dim;
	struct shape *s;
	long double x_pos=0.0;
	long double y_pos=0.0;
	long double z_pos=0.0;

	for (z=0;z<dim->zlen;z++)
	{
		z_pos=dim->zmesh[z];
		for (x=0;x<dim->xlen;x++)
		{
			x_pos=dim->xmesh[x];
			for (y=0;y<dim->ylen;y++)
			{
				y_pos=dim->ymesh[y];
				for (c=0;c<in->ncontacts;c++)
				{
					if (in->contacts[c].position==LEFT)		//I should not need this line if I fix up the other contacts
					{
						s=&(in->contacts[c].shape);
						if (shape_in_shape(sim,s,z_pos,x_pos,y_pos+in->my_epitaxy.device_start)==0)
						{
							in->mun[z][x][y]=1.0;
							in->mup[z][x][y]=1.0;
						}
					}

				}

				//printf("%d %d %Le\n",x,y,in->mun[z][x][y]);
			}
		}

	}

//getchar();
}
