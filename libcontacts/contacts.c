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
#include "const.h"
#include <cal_path.h>
#include "contacts.h"
#include <log.h>
#include <dump.h>

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

			if (in->last_ittr<16)
			{
				dV=dV*(1.0+0.001);
			}

			if (in->last_ittr>18)
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

			printf_log(sim,"Ramping: %s %.2Lf %.2Lf dV=%Lf ittr=%d\n",in->contacts[i].name,in->contacts[i].voltage,in->contacts[i].voltage_want,dV,in->last_ittr);

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

	inp_check(sim,&inp,1.2);
	inp_reset_read(sim,&inp);
	inp_get_string(sim,&inp);
	sscanf(inp_get_string(sim,&inp),"%d",&(in->ncontacts));

	if (in->ncontacts>10)
	{
		ewe(sim,"Too many contacts\n");
	}

	if (in->ncontacts<1)
	{
		ewe(sim,"No contacts\n");
	}

	gdouble pos=0.0;
	int active=FALSE;
	for (i=0;i<in->ncontacts;i++)
	{

		inp_get_string(sim,&inp);	//active contact
		strcpy(in->contacts[i].name,inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//position
		in->contacts[i].position=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//active contact
		in->contacts[i].active=english_to_bin(sim, inp_get_string(sim,&inp));

		inp_get_string(sim,&inp);	//start
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].start));

		inp_get_string(sim,&inp);	//width
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].width));

		inp_get_string(sim,&inp);	//depth
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].depth));

		inp_get_string(sim,&inp);	//voltage
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].voltage_want));
		in->contacts[i].voltage=0.0;
		in->contacts[i].voltage_last=in->contacts[i].voltage;

		inp_get_string(sim,&inp);	//np
		sscanf(inp_get_string(sim,&inp),"%Le",&(in->contacts[i].np));
		in->contacts[i].np=fabs(in->contacts[i].np);

		inp_get_string(sim,&inp);	//charge_type
		in->contacts[i].charge_type=english_to_bin(sim, inp_get_string(sim,&inp));

		pos+=in->contacts[i].width;
	}

	char * ver = inp_get_string(sim,&inp);
	if (strcmp(ver,"#ver")!=0)
	{
			ewe(sim,"No #ver tag found in file\n");
	}

	inp_free(sim,&inp);

	contacts_update(sim,in);
	contact_set_flip_current(sim,in);

	in->lcharge=contacts_get_lcharge(sim,in);
	in->rcharge=contacts_get_rcharge(sim,in);


	contacts_cal_area(sim,in);
	//contacts_dump(sim,in);
	//getchar();
}

void contacts_force_to_zero(struct simulation *sim,struct device *in)
{
int x;
int z;

for (x=0;x<in->xmeshpoints;x++)
{
	for (z=0;z<in->zmeshpoints;z++)
	{
		in->Vapplied_l[z][x]=0.0;
		in->Vapplied_r[z][x]=0.0;
	}

}

}


void contacts_dump(struct simulation *sim,struct device *in)
{
int i;
	if (get_dump_status(sim,dump_print_text)==TRUE)
	{
		for (i=0;i<in->ncontacts;i++)
		{
			printf_log(sim,"%-10s\tV=%Le\tA=%Le\n",in->contacts[i].name,in->contacts[i].voltage,in->contacts[i].area);
		}


		//printf("%Le\n",in->flip_current);
	}

	int z=0;
	int x=0;


	for (z = 0; z < in->zmeshpoints; z++)
	{
		for (x = 0; x < in->xmeshpoints; x++)
		{
			printf("%.2Lf %.2Lf\n",in->Vapplied_l[z][x],in->Vapplied_r[z][x]);
		}
	}

}

void contacts_update(struct simulation *sim,struct device *in)
{
int i;
int x;
int z;
int n;
int found=FALSE;

gdouble value=0.0;

if (in->xmeshpoints==1)
{
	for (z=0;z<in->zmeshpoints;z++)
	{
		in->passivate_r[z][0]=FALSE;
		in->passivate_l[z][0]=FALSE;
		for (i=0;i<in->ncontacts;i++)
		{
			if ((in->contacts[i].position==TOP)&&(in->contacts[i].active==TRUE))
			{
				in->Vapplied_l[z][0]=in->contacts[i].voltage;
			}else
			{
				in->Vapplied_r[z][0]=in->contacts[i].voltage;
			}
		}
	}

	//contacts_dump(sim,in);
	//getchar();
	return;
}

//Reset contacts
for (z=0;z<in->zmeshpoints;z++)
{
	for (x=0;x<in->xmeshpoints;x++)
	{
		in->Vapplied_l[z][x]=0.0;
		in->Vapplied_r[z][x]=0.0;

		in->n_contact_r[z][x]=-1;
		in->n_contact_l[z][x]=-1;

		in->passivate_r[z][x]=TRUE;
		in->passivate_l[z][x]=TRUE;
	}
}

for (i=0;i<in->ncontacts;i++)
{
	for (z=0;z<in->zmeshpoints;z++)
	{
		for (x=0;x<in->xmeshpoints;x++)
		{
			if ((in->xmesh[x]>=in->contacts[i].start)&&(in->xmesh[x]<in->contacts[i].start+in->contacts[i].width))
			{

				if (in->contacts[i].position==TOP)
				{
					if (in->n_contact_l[z][x]!=-1)
					{
						ewe(sim,"You have overlapping contacts\n");
					}

					in->Vapplied_l[z][x]=in->contacts[i].voltage;
					in->n_contact_l[z][x]=i;
					in->passivate_l[z][x]=FALSE;
				}else
				{
					if (in->n_contact_r[z][x]!=-1)
					{
						ewe(sim,"You have overlapping contacts\n");
					}

					in->Vapplied_r[z][x]=in->contacts[i].voltage;
					in->n_contact_r[z][x]=i;
					in->passivate_r[z][x]=FALSE;
				}
			}
		}
	}


}

//contacts_dump(sim,in);

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

for (x=0;x<in->xmeshpoints;x++)
{
		for (z=0;z<in->zmeshpoints;z++)
		{
			if (in->n_contact_l[z][x]>=0)
			{
				tot+=in->Jpleft[z][x]+in->Jnleft[z][x];
				count=count+1.0;						//this will need updating for meshes which change
			}
		}
}

tot=tot/count;
tot*=in->flip_current;

return tot*Q;
}

long double contacts_get_Jright(struct device *in)
{
int i;
int x;
int z;

long double tot=0.0;
long double count=0.0;

for (x=0;x<in->xmeshpoints;x++)
{
		for (z=0;z<in->zmeshpoints;z++)
		{
			if (in->n_contact_r[z][x]>=0)
			{
				tot+=in->Jpright[z][x]+in->Jnright[z][x];
				count=count+1.0;
			}
		}
}

tot=tot/count;

tot*=in->flip_current;

return tot*Q;
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

for (i=0;i<in->ncontacts;i++)
{
	in->contacts[i].area=0.0;
}

for (x=0;x<in->xmeshpoints;x++)
{
		for (z=0;z<in->zmeshpoints;z++)
		{
			i=in->n_contact_r[z][x];
			if (i!=-1)
			{
				in->contacts[i].area+=in->dxmesh[x]*in->dzmesh[z];
			}

			i=in->n_contact_l[z][x];
			if (i!=-1)
			{
				in->contacts[i].area+=in->dxmesh[x]*in->dzmesh[z];
			}

		}
}

}

void contacts_detailed_dump(struct device *in)
{
int i;
int x;
int z;

for (x=0;x<in->xmeshpoints;x++)
{
		for (z=0;z<in->zmeshpoints;z++)
		{
			printf("%d %Le %Le\n",in->n_contact_l[z][x],in->Jpleft[z][x],in->Jnleft[z][x]);
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

for (x=0;x<in->xmeshpoints;x++)
{
		for (z=0;z<in->zmeshpoints;z++)
		{
			for (i=0;i<in->ncontacts;i++)
			{
				if (in->n_contact_r[z][x]==n)
				{
					tot+=in->Jpright[z][x]+in->Jnright[z][x];
					count=count+1.0;						//this will need updating for meshes which change
				}

				if (in->n_contact_l[z][x]==n)
				{
					tot+=in->Jpleft[z][x]+in->Jnleft[z][x];
					count=count+1.0;						//this will need updating for meshes which change
				}
			}
		}
}

tot=tot/count;

tot*=in->flip_current;

return tot*Q;
}

void contacts_passivate(struct simulation *sim,struct device *in)
{
int i;
int x;
int y;
int z;
return;
//passivate under each contact
for (x=0;x<in->xmeshpoints;x++)
{
	for (y=0;y<in->ymeshpoints;y++)
	{
		for (z=0;z<in->zmeshpoints;z++)
		{

			for (i=0;i<in->ncontacts;i++)
			{
				if (in->contacts[i].position==TOP)
				{
					if ((in->ylen-in->ymesh[y]<=in->contacts[i].depth)&&(in->xmesh[x]>in->contacts[i].start)&&(in->xmesh[x]<in->contacts[i].start+in->contacts[i].width))
					{
						in->mun[z][x][y]=1e-15;
						in->mup[z][x][y]=1e-15;
					}
				}

				if (in->contacts[i].position==BOTTOM)
				{
					if ((in->ymesh[y]<=in->contacts[i].depth)&&(in->xmesh[x]>in->contacts[i].start)&&(in->xmesh[x]<in->contacts[i].start+in->contacts[i].width))
					{
						in->mun[z][x][y]=1e-15;
						in->mup[z][x][y]=1e-15;
					}
				}
			}
		}
	}
}

for (x=0;x<in->xmeshpoints;x++)
{
	for (z=0;z<in->zmeshpoints;z++)
	{
		i=in->n_contact_r[z][x];
		if (i==-1)
		{
			in->mun[z][x][in->ymeshpoints-1]=1e-15;
			in->mup[z][x][in->ymeshpoints-1]=1e-15;
		}

		i=in->n_contact_l[z][x];
		if (i==-1)
		{
			in->mun[z][x][0]=1e-15;
			in->mup[z][x][0]=1e-15;
		}
	}
}

}
